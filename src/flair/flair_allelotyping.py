#! /usr/bin/env python3

import argparse
import pysam
from flair.gtf_io import load_gtf_to_gene_data, GtfExon
from copy import deepcopy
import graphviz
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
import math
import vcfpy
import networkx as nx
from collections import OrderedDict
import string

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bam', required=True, type=str,
                        help="bam file (this is your primary bam input, required)")
    parser.add_argument('--norm_bam', type=str,
                        help="normal bam file (optional, used to call variants as somatic or not)")
    parser.add_argument('--vcf', required=True, type=str,
                        help="vcf file (even if you are inputting both normal and tumor bams, you can still provide a single vcf here)")
    parser.add_argument('--norm_vcf', type=str,
                        help="normal vcf file, optional even if providing a normal bam")
    parser.add_argument('-o', '--output', default='flair',
                        help="output prefix. default: 'flair'")
    parser.add_argument('-g', '--genome',
                        type=str, required=True, help='FastA of reference genome')
    parser.add_argument('-f', '--gtf', type=str,
                        help='GTF annotation file, either this or isoform_bed must be provided')
    parser.add_argument('--isoform_bed', type=str,
                        help="FLAIR isoform bed file, can be used as primary or supplementary annotation")
    parser.add_argument('--frac_support', default=0.05, type=float,
                        help='fraction of reads required to call allele, not recommended to change from default')
    parser.add_argument('--read_support', default=3, type=int,
                        help='number of reads required to call allele')
    parser.add_argument('--annotate_bams', default=False, action='store_true',
                        help='label reads in bam files with HP tags for allele groups')
    parser.add_argument('--generate_map', default=False, action='store_true',
                        help='output read map of allele groups to read names')
    # parser.add_argument('--keep_intermediate', default=False, action='store_true',
    #                     help='''specify if intermediate and temporary files are to be kept for debugging.''')
    args = parser.parse_args()
    if args.gtf is None and args.isoform_bed is None:
        parser.error('At least one of gtf and isoform_bed must be provided')
    return args

def make_vizgraph(graph, readvarinfo_to_reads):
    vizgraph = graphviz.Digraph()
    c = 0
    og_node_to_index = {}
    for node in graph.nodes:
        strnode = [f"{x[0]};{x[1]}" for x in node]
        vizgraph.node(str(c), f'{", ".join(strnode)}: {len(readvarinfo_to_reads[node])}')
        og_node_to_index[node] = str(c)
        c += 1
    for edge in graph.edges:
        vizgraph.edge(og_node_to_index[edge[0]], og_node_to_index[edge[1]])
    vizgraph.render('ULC0515T-variants-graph-simpler_KRAS_060426.gv')

def combine_vcf_files(vcffilelist):
    vartoalt = {}
    for samplevcf in vcffilelist:
        for line in open(samplevcf):
            if line[0] != '#':
                line = line.rstrip().split('\t')
                if line[6] == 'PASS':
                    ref = line[3]
                    refinfo, alts = (line[0], int(line[1]) - 1), line[4].split(',')
                    if refinfo not in vartoalt:
                        vartoalt[refinfo] = set()
                    for alt in alts:
                        vartoalt[refinfo].add((ref, alt))
    return vartoalt

def reorganize_vars(vartoalt):
    chrom_to_region_to_vars = {}
    for chrom, pos in vartoalt:
        vars = tuple(vartoalt[(chrom, pos)])
        region = math.floor(pos / (10 ** 5))
        if chrom not in chrom_to_region_to_vars:
            chrom_to_region_to_vars[chrom] = {}
        if region not in chrom_to_region_to_vars[chrom]:
            chrom_to_region_to_vars[chrom][region] = set()
        chrom_to_region_to_vars[chrom][region].add((pos, vars))
    return chrom_to_region_to_vars

def get_gene_to_all_vars(gene_to_all_exons, chrom_to_region_to_vars, gene_to_chrom=None):
    gene_to_vars = {}
    for gene in gene_to_all_exons:
        if gene_to_chrom is not None:
            chrom = gene_to_chrom[gene]
        else:
            chrom = gene_to_all_exons[gene][0].chrom
        gene_to_vars[gene] = {}
        if chrom in chrom_to_region_to_vars:
            for block in gene_to_all_exons[gene]:
                for region in range(math.floor(block.start / (10 ** 5)), math.floor(block.end / (10 ** 5)) + 1):
                    if region in chrom_to_region_to_vars[chrom]:
                        for pos, vars in chrom_to_region_to_vars[chrom][region]:
                            if block.start <= pos <= block.end:
                                gene_to_vars[gene][pos] = vars
    return gene_to_vars

def get_gene_boundaries(exons, chrom=None):
    if chrom is None:
        chrom = exons[0].chrom
    start = min([x.start for x in exons])
    end = max([x.end for x in exons])
    return chrom, start, end

def parse_cigar_for_introns_indels(cigar, alignstart):
    introns, insertions, deletions = [], [], []
    ref, quer = alignstart, 0
    for block in cigar:
        if block[0] in {0, 7, 8}:  # match, consumes both
            ref += block[1]
            quer += block[1]
        elif block[0] in {1, 4}:  # consumes query, 1 is insertion
            if block[0] == pysam.CIGAR_OPS.CINS:
                insertions.append((ref, quer, block[1]))
            quer += block[1]
        elif block[0] in {2, 3}:  # consumes reference ##2 is deletion
            if block[0] == pysam.CIGAR_OPS.CDEL:
                deletions.append((ref, quer, block[1]))
            elif block[0] == pysam.CIGAR_OPS.CREF_SKIP:
                introns.append((ref, ref + block[1]))
            ref += block[1]
    return introns, insertions, deletions

def add_indels_to_vars(insertions, deletions, chrom, readseq, genome, my_vars):
    for refpos, querpos, length in insertions:
        if refpos - 1 in my_vars:
            vinfo = (readseq[querpos - 1], readseq[querpos - 1:querpos + length])
            if vinfo in my_vars[refpos - 1]:
                my_vars[refpos - 1][vinfo][1] = 1
    for refpos, querpos, length in deletions:
        if refpos - 1 in my_vars:
            vinfo = (genome.fetch(chrom, refpos - 1, refpos + length).upper(), genome.fetch(chrom, refpos - 1, refpos).upper())
            if vinfo in my_vars[refpos - 1]:
                my_vars[refpos - 1][vinfo][1] = 1

def get_coverage_snvs_from_read(a, my_vars):
    alignedbases = a.get_aligned_pairs(with_seq=True, matches_only=True)
    for querpos, refpos, base in alignedbases:
        if refpos in my_vars:
            for v in my_vars[refpos]:
                my_vars[refpos][v][0] = 1
            if base.islower():
                vinfo = (base.upper(), a.query_sequence[querpos].upper())
                if vinfo in my_vars[refpos]:
                    my_vars[refpos][vinfo][1] = 1

def check_any_var_is_covered(my_vars):
    for pos in my_vars:
        for ref, alt in my_vars[pos]:
            cov, var = my_vars[pos][ref, alt]
            if cov != 0:
                return True
    return False

def load_bam_read_vars(a, bam_label, my_vars, genome, insertions, deletions):
    add_indels_to_vars(insertions, deletions, a.reference_name, a.query_sequence, genome, my_vars)
    get_coverage_snvs_from_read(a, my_vars)
    if check_any_var_is_covered(my_vars):
        posinfo, thisreadinfo = simplify_gene_vars(my_vars)
        return (bam_label, a.query_name), thisreadinfo
    return None, None

def create_gene_vars_dict(gene_vars):
    my_vars = {}
    for pos in gene_vars:
        my_vars[pos] = {}
        for (ref, alt) in gene_vars[pos]:
            my_vars[pos][(ref, alt)] = [0, 0]
    return my_vars

def simplify_gene_vars(my_vars):
    info = []
    for pos in my_vars:
        for ref, alt in my_vars[pos]:
            cov, var = my_vars[pos][ref, alt]
            info.append((pos, ref, alt, cov, var))
    info.sort()
    posinfo = [f'{x[0]};{x[1]};{x[2]}' for x in info]
    thisreadinfo = []
    for index, (pos, ref, alt, cov, var) in enumerate(info):
        if cov == 1:
            thisreadinfo.append((index, var))
    return posinfo, tuple(thisreadinfo)

def add_nodes(graph, readvarinfo_to_reads):
    tot_gene_reads = 0
    for a in readvarinfo_to_reads:
        graph.add_node(a)
        tot_gene_reads += len(readvarinfo_to_reads[a])
    return tot_gene_reads

def add_edges_from_var_supersets(graph, readvarinfo_to_reads):
    # construct edges based on variant supersets
    for a in readvarinfo_to_reads:
        for b in readvarinfo_to_reads:
            if a != b and len(set(a) - set(b)) == 0:
                graph.add_edge(a, b)

def simplify_remove_nodes_with_single_parent(graph, readvarinfo_to_reads):
    g2 = graph.copy()
    moved_nodes = {}
    for node in graph:
        if len(list(graph.out_edges(node))) == 1:
            thisnode, nextnode = list(graph.out_edges(node))[0]
            while nextnode in moved_nodes:
                nextnode = moved_nodes[nextnode]
            # update edges upstream
            for prev, this in list(graph.in_edges(node)):
                g2.add_edge(prev, nextnode)
            readvarinfo_to_reads[nextnode].update(readvarinfo_to_reads[thisnode])
            readvarinfo_to_reads[thisnode] = set()
            moved_nodes[thisnode] = nextnode
            g2.remove_node(node)
    return g2

def remove_low_support_terminal_nodes(graph, readvarinfo_to_reads, tot_gene_reads, reads_lost, do_frac, read_support, frac_support):
    # remove terminal nodes that don't meet fraction or coverage standards
    g2 = graph.copy()
    for node in graph:
        if len(list(graph.out_edges(node))) == 0:
            if (do_frac and len(readvarinfo_to_reads[node]) / tot_gene_reads < frac_support) \
                    or (not do_frac and len(readvarinfo_to_reads[node]) < read_support):
                reads_lost.update(readvarinfo_to_reads[node])
                readvarinfo_to_reads[node] = set()
                g2.remove_node(node)
    return g2

def get_terminal_node(graph, node, terminal_nodes):
    if len(list(graph.out_edges(node))) > 0:
        for thisnode, nextnode in graph.out_edges(node):
            terminal_nodes.update(get_terminal_node(graph, nextnode, terminal_nodes))
    else:
        terminal_nodes.add(node)
    return terminal_nodes

def assign_support_to_terminal_nodes(graph, readvarinfo_to_reads, reads_lost, remove_ambig):
    g2 = graph.copy()
    for node in graph:
        if len(list(graph.out_edges(node))) != 0:
            term_nodes = list(get_terminal_node(graph, node, set()))
            if len(term_nodes) == 1:
                readvarinfo_to_reads[term_nodes[0]].update(readvarinfo_to_reads[node])
                g2.remove_node(node)
                readvarinfo_to_reads[node] = set()
            elif remove_ambig:
                reads_lost.update(readvarinfo_to_reads[node])
                readvarinfo_to_reads[node] = set()
                g2.remove_node(node)
    return g2

def identify_final_allele_groups(graph, readvarinfo_to_reads):
    allele_group_to_reads = {}
    for node in graph.nodes:
        allele_group_to_reads[node] = readvarinfo_to_reads[node]
    return allele_group_to_reads


def process_alleotype_graph(readvarinfo_to_reads, read_support, frac_support):
    reads_lost = set()
    graph = nx.DiGraph()
    tot_gene_reads = add_nodes(graph, readvarinfo_to_reads)
    add_edges_from_var_supersets(graph, readvarinfo_to_reads)

    graph = simplify_remove_nodes_with_single_parent(graph, readvarinfo_to_reads)

    graph = remove_low_support_terminal_nodes(graph, readvarinfo_to_reads, tot_gene_reads, reads_lost, False, read_support, frac_support)
    graph = assign_support_to_terminal_nodes(graph, readvarinfo_to_reads, reads_lost, remove_ambig=False)

    graph = remove_low_support_terminal_nodes(graph, readvarinfo_to_reads, tot_gene_reads, reads_lost, True, read_support, frac_support)
    graph = assign_support_to_terminal_nodes(graph, readvarinfo_to_reads, reads_lost, remove_ambig=True)

    return identify_final_allele_groups(graph, readvarinfo_to_reads)

def label_bam_file(bamname, output_name, read_to_allele_group):
    with pysam.AlignmentFile(bamname, 'rb') as bam:
        with pysam.AlignmentFile(output_name, 'wb', template=bam) as out:
            for a in bam:
                if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                    if a.query_name in read_to_allele_group:
                        a.set_tag('HP', read_to_allele_group[a.query_name])
                out.write(a)
    pysam.index(output_name)

def label_bam_files(bam, norm_bam, output, file_to_read_to_allele_group):
    label_bam_file(bam, output + '.tumor.bam', file_to_read_to_allele_group['t'])
    if norm_bam is not None:
        label_bam_file(norm_bam, output + '.normal.bam', file_to_read_to_allele_group['n'])

def generate_new_vcf_header(in_vcf, norm_bam):
    with vcfpy.Reader.from_path(in_vcf) as reader:
        new_header = vcfpy.Header()
        for line in reader.header.lines:
            if type(line) in (vcfpy.HeaderLine, vcfpy.ContigHeaderLine):
                new_header.add_line(line)
    new_header.add_line(vcfpy.HeaderLine('source', 'FLAIR allelotyping'))
    new_header.add_line(vcfpy.InfoHeaderLine.from_mapping({'ID': 'PS', 'Number': 1, 'Type': 'Integer', 'Description': 'Phase set - variants within a phase set can be compared/phased'}))
    new_header.add_line(vcfpy.InfoHeaderLine.from_mapping({'ID': 'AG', 'Number': 'R', 'Type': 'String', 'Description': 'Allele group(s) - only compare these within a phase set'}))
    new_header.add_line(vcfpy.InfoHeaderLine.from_mapping({'ID': 'DF', 'Number': 1, 'Type': 'Integer', 'Description': '1 means this is a variant that helps differentiate this phase set from others'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'DP', 'Number': 1, 'Type': 'Integer', 'Description': 'Read depth'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths for the ref and alt alleles in the order listed'}))
    if norm_bam is not None:
        new_header.samples = vcfpy.SamplesInfos(['tumor', 'normal'])
    else:
        new_header.samples = vcfpy.SamplesInfos(['tumor'])
    return new_header

def write_vcf_file(output, new_header, variant_to_allele_group_counts_info, norm_bam, phaseset_to_allele_groups):
    allele_group_to_final_vars = {}
    for ps in phaseset_to_allele_groups:
        for ag in phaseset_to_allele_groups[ps]:
            allele_group_to_final_vars[(ps, ag)] = set()
    with vcfpy.Writer.from_path(output + '.allelegroups.vcf', new_header) as writer:
        format_strings = ['GT', 'DP', 'AD']
        for (chrom, pos, ref, alt), var_data in variant_to_allele_group_counts_info.items():
            tot_cov = var_data['t_cov'] + var_data['n_cov']
            tot_var = var_data['t_var'] + var_data['n_var']
            # NEW: reporting homozygous variants + variants in all allele groups - more important for downstream vars
            if tot_var != 0 and len(var_data['ag']) != 0:
                is_defining_var = tot_var != tot_cov and len(var_data['ag']) != len(phaseset_to_allele_groups[var_data['ag'][0][0]])
                gt = '0/1' if is_defining_var else '1/1'
                is_defining_var = 1 if is_defining_var else 0
                # FIXME allow other types of alts
                alt_type = 'SNV'
                alt_desc = vcfpy.Substitution(alt_type, alt)
                these_calls = [vcfpy.Call('tumor', {'GT': gt, 'DP': var_data['t_cov'], 'AD': [var_data['t_cov'] - var_data['t_var'], var_data['t_var']]})]
                if norm_bam is not None:
                    these_calls.append(vcfpy.Call('normal', {'GT': gt, 'DP': var_data['n_cov'], 'AD': [var_data['n_cov'] - var_data['n_var'], var_data['n_var']]}))
                my_ps = list(set([x[0] for x in var_data['ag']]))
                if len(my_ps) > 1:
                    raise ValueError(f'multiple phase sets for one variant {my_ps} {var_data}')
                my_ag = sorted([x[1] for x in var_data['ag']])
                new_record = vcfpy.Record(chrom, int(pos) + 1, [], ref, [alt_desc], 100, ['PASS'], OrderedDict([('PS', my_ps[0]), ('AG', my_ag), ('DF', is_defining_var)]), format_strings, these_calls)
                writer.write_record(new_record)
                for allele_group in var_data['ag']:
                    # if allele_group not in allele_group_to_final_vars:
                    #     allele_group_to_final_vars[allele_group] = set()
                    allele_group_to_final_vars[allele_group].add((chrom, int(pos), ref, alt))
    return allele_group_to_final_vars

def write_allele_group_counts_read_map(index_to_allele_group_info, output, generate_map, norm_bam):
    with open(output + '.allelegroups.counts.tsv', 'w') as fh:
        outline = ['#phase_set', 'allele_group', 'gene', 'tumor_counts']
        if norm_bam is not None:
            outline.append('normal_counts')
        fh.write('\t'.join(outline) + '\n')
        for (ps, ag), group_info in index_to_allele_group_info.items():
            ol = [str(ps), ag, group_info['gene'], str(len(group_info['t']))]
            if norm_bam is not None:
                ol.append(str(len(group_info['n'])))
            fh.write('\t'.join(ol) + '\n')
    if generate_map:
        with open(output + '.allelegroups.tumor.read.map.txt', 'w') as fh:
            for (ps, ag), group_info in index_to_allele_group_info.items():
                fh.write(f'{ps}|{ag}\t{",".join(sorted(group_info["t"]))}\n')
        if norm_bam is not None:
            with open(output + '.allelegroups.normal.read.map.txt', 'w') as fh:
                for (ps, ag), group_info in index_to_allele_group_info.items():
                    fh.write(f'{ps}|{ag}\t{",".join(sorted(group_info["n"]))}\n')

def make_allele_group_label(phaseset, group, allele_group_count, allele_group_to_reads, norm_bam, read_support):
    tot_reads_for_file = {'t': 0, 'n': 0}
    for file_label, read in allele_group_to_reads[group]:
        tot_reads_for_file[file_label] += 1
    allele_group_label = (phaseset, string.ascii_uppercase[allele_group_count])
    # identify if is somatic
    if norm_bam is not None and tot_reads_for_file['n'] < read_support:
        allele_group_label = (allele_group_label[0], allele_group_label[1] + '-S')
    allele_group_count += 1
    return allele_group_count, allele_group_label

def get_phase_sets(allele_group_to_reads, var_info, phaseset_count):
    pos_ranges = []
    for group in allele_group_to_reads:
        covered_pos = []
        for varindex, is_var in group:
            covered_pos.append(int(var_info[varindex].split(';')[0]))
        pos_ranges.append((min(covered_pos), max(covered_pos), group))
    pos_ranges.sort()
    phaseset_to_groups = {}
    lastend, grouped = -1, []
    for start, end, group in pos_ranges:
        if start > lastend:
            if len(grouped) > 0:
                phaseset_to_groups[phaseset_count] = grouped
                phaseset_count += 1
                grouped = []
        if start > lastend or end > lastend:
            lastend = end
        grouped.append(group)
    if len(grouped) > 0:
        phaseset_to_groups[phaseset_count] = grouped
        phaseset_count += 1
    return phaseset_to_groups, phaseset_count

def get_allele_group_info(allele_group_to_reads, phaseset_count, var_info, gene_id, gene_chrom, index_to_allele_group_info, file_to_read_to_allele_group, variant_to_allele_group_counts_info, phaseset_to_allele_groups, norm_bam, read_support):
    phaseset_to_groups, phaseset_count = get_phase_sets(allele_group_to_reads, var_info, phaseset_count)
    for ps in phaseset_to_groups:
        allele_group_count = 0
        phaseset_to_allele_groups[ps] = set()
        for group in phaseset_to_groups[ps]:
            allele_group_count, allele_group_label = make_allele_group_label(ps, group, allele_group_count, allele_group_to_reads, norm_bam, read_support)
            phaseset_to_allele_groups[ps].add(allele_group_label[1])
            index_to_allele_group_info[allele_group_label] = {'gene': gene_id, 'vars': var_info, 'has_vars': group, 't': set(), 'n': set()}
            for file_label, read in allele_group_to_reads[group]:
                file_to_read_to_allele_group[file_label][read] = '|'.join([str(x) for x in allele_group_label])
                index_to_allele_group_info[allele_group_label][file_label].add(read)

            for varindex, is_var in group:
                pos, ref, alt = var_info[varindex].split(';')
                var_data = (gene_chrom, pos, ref, alt)
                if var_data not in variant_to_allele_group_counts_info:
                    variant_to_allele_group_counts_info[var_data] = {'gene': gene_id, 'ag': [], 't_cov': 0, 't_var': 0, 'n_cov': 0, 'n_var': 0}
                if is_var == 1:
                    variant_to_allele_group_counts_info[var_data]['ag'].append(allele_group_label)

    return phaseset_count

def get_variant_final_coverage(readvarinfo_to_reads, var_info, gene_chrom, variant_to_allele_group_counts_info):
    for og_group in readvarinfo_to_reads:
        for varindex, is_var in og_group:
            pos, ref, alt = var_info[varindex].split(';')
            var_data = (gene_chrom, pos, ref, alt)
            if var_data in variant_to_allele_group_counts_info:
                for file_label, read in readvarinfo_to_reads[og_group]:
                    variant_to_allele_group_counts_info[var_data][file_label + '_cov'] += 1
                    if is_var == 1:
                        variant_to_allele_group_counts_info[var_data][file_label + '_var'] += 1

def check_read_is_in_gene(gene_juncs, introns):
    if len(set(gene_juncs) & set(introns)) > 0:
        return True
    return False

def load_bam_file_for_region(bam, bam_label, rchrom, rstart, rend, my_vars, genome, gene_juncs, readvarinfo_to_reads):
    with pysam.AlignmentFile(bam, 'rb') as bam:
        for a in bam.fetch(rchrom, rstart, rend):
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                introns, insertions, deletions = parse_cigar_for_introns_indels(a.cigartuples, a.reference_start)
                if check_read_is_in_gene(gene_juncs, introns):
                    read_id, read_vars = load_bam_read_vars(a, bam_label, deepcopy(my_vars), genome, insertions, deletions)
                    if read_vars is not None:
                        if read_vars not in readvarinfo_to_reads:
                            readvarinfo_to_reads[read_vars] = set()
                        readvarinfo_to_reads[read_vars].add(read_id)

def load_bam_files_for_region(bam, norm_bam, gene_chrom, gene_start, gene_end, my_vars, genome, gene_juncs):
    readvarinfo_to_reads = {}
    load_bam_file_for_region(bam, 't', gene_chrom, gene_start, gene_end, my_vars, genome, gene_juncs, readvarinfo_to_reads)
    if norm_bam is not None:
        load_bam_file_for_region(norm_bam, 'n', gene_chrom, gene_start, gene_end, my_vars, genome, gene_juncs, readvarinfo_to_reads)
    return readvarinfo_to_reads

def load_isoform_bed_data(gtf, isoform_bed, gene_to_all_exons, gene_to_all_juncs):
    for bed in BedReader(isoform_bed, bedClass=FlairBed):
        if bed.transcript_class != 'fusion':
            gene = None
            if gtf is None or len(bed.ref_gene_mappings) == 0:
                gene = bed.gene_id
            elif len(bed.ref_gene_mappings) == 1:
                gene = bed.ref_gene_mappings[0]
            if gene is not None:
                if gene not in gene_to_all_exons:
                    gene_to_all_exons[gene] = set()
                    gene_to_all_juncs[gene] = set()
                for exon in bed.blocks:
                    gene_to_all_exons[gene].add(GtfExon(bed.chrom, None, 'exon', exon.start, exon.end, None, bed.strand, None))
                for i in range(len(bed.blocks) - 1):
                    gene_to_all_juncs[gene].add((bed.blocks[i].end, bed.blocks[i + 1].start))

def load_isoform_data(gtf, isoform_bed):
    gene_to_all_exons, gene_to_all_juncs = {}, {}
    if gtf is not None:
        gene_to_all_exons, gene_to_all_juncs = load_gtf_to_gene_data(gtf)
    if isoform_bed is not None:
        load_isoform_bed_data(gtf, isoform_bed, gene_to_all_exons, gene_to_all_juncs)
    for gene in gene_to_all_exons:
        gene_to_all_exons[gene] = sorted(list(gene_to_all_exons[gene]), key=lambda x: (x.start, x.end))
        gene_to_all_juncs[gene] = sorted(list(gene_to_all_juncs[gene]))
    return gene_to_all_exons, gene_to_all_juncs

def getvariants():
    args = parse_args()
    print('loading genes and variants')
    genome = pysam.FastaFile(args.genome)
    my_vcfs = [args.vcf, args.norm_vcf] if args.norm_vcf is not None else [args.vcf]
    vartoalt = combine_vcf_files(my_vcfs)
    chrom_to_region_to_vars = reorganize_vars(vartoalt)
    gene_to_all_exons, gene_to_all_juncs = load_isoform_data(args.gtf, args.isoform_bed)

    gene_to_vars = get_gene_to_all_vars(gene_to_all_exons, chrom_to_region_to_vars)
    # tempdir = make_temp_dir(args.output)

    print('loading bam files and creating allele groups')
    file_to_read_to_allele_group = {'t': {}, 'n': {}}
    index_to_allele_group_info = {}
    variant_to_allele_group_counts_info = {}
    phaseset_to_allele_groups = {}
    phaseset_count = 1
    # TODO can parallelize at the gene level, will need to write to intermediate files
    for gene_id in gene_to_vars:
        # only process reads if gene has variants
        # and gene is spliced (junctions are how we assign reads to gene)
        if len(gene_to_vars[gene_id]) > 0 and len(gene_to_all_juncs[gene_id]) > 0:
            gene_chrom, gene_start, gene_end = get_gene_boundaries(gene_to_all_exons[gene_id])
            my_vars = create_gene_vars_dict(gene_to_vars[gene_id])
            var_info, _ = simplify_gene_vars(my_vars)
            readvarinfo_to_reads = load_bam_files_for_region(args.bam, args.norm_bam, gene_chrom, gene_start, gene_end, my_vars, genome, gene_to_all_juncs[gene_id])
            allele_group_to_reads = process_alleotype_graph(deepcopy(readvarinfo_to_reads), args.read_support, args.frac_support)

            # don't report if there's only one final group
            # if len(allele_group_to_reads) > 1: allele_group_to_reads, phaseset_count, var_info, gene_id
            phaseset_count = get_allele_group_info(allele_group_to_reads, phaseset_count, var_info, gene_id, gene_chrom, index_to_allele_group_info,
                                                   file_to_read_to_allele_group, variant_to_allele_group_counts_info, phaseset_to_allele_groups, args.norm_bam, args.read_support)

            # getting variant coverage from original assignments, not collapsed ones
            get_variant_final_coverage(readvarinfo_to_reads, var_info, gene_chrom, variant_to_allele_group_counts_info)

    print('writing vcf file')
    new_header = generate_new_vcf_header(args.vcf, args.norm_bam)
    write_vcf_file(args.output, new_header, variant_to_allele_group_counts_info, args.norm_bam, phaseset_to_allele_groups)

    write_allele_group_counts_read_map(index_to_allele_group_info, args.output, args.generate_map, args.norm_bam)

    if args.annotate_bams:
        print('labeling bam files')
        label_bam_files(args.bam, args.norm_bam, args.output + '.allelegroups', file_to_read_to_allele_group)


if __name__ == "__main__":
    getvariants()
