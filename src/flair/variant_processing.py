#! /usr/bin/env python3

import math
import vcfpy
import pysam
import networkx as nx
from collections import OrderedDict

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
        if refpos in my_vars:
            vinfo = (readseq[querpos - 1], readseq[querpos - 1:querpos + length])
            if vinfo in my_vars:
                my_vars[refpos][vinfo][1] = 1
    for refpos, querpos, length in deletions:
        if refpos in my_vars:
            vinfo = (refpos, genome.fetch(chrom, refpos - 1, refpos + length).upper(), genome.fetch(chrom, refpos - 1, refpos + length).upper())
            if vinfo in my_vars:
                my_vars[refpos][vinfo][1] = 1

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

def load_bam_read_vars(a, bam_label, my_vars, genome, insertions, deletions):
    add_indels_to_vars(insertions, deletions, a.reference_name, a.query_sequence, genome, my_vars)
    get_coverage_snvs_from_read(a, my_vars)
    if check_any_var_is_covered(my_vars):
        posinfo, thisreadinfo = simplify_gene_vars(my_vars)
        return (bam_label, a.query_name), thisreadinfo
    return None, None

def check_any_var_is_covered(my_vars):
    for pos in my_vars:
        for ref, alt in my_vars[pos]:
            cov, var = my_vars[pos][ref, alt]
            if cov != 0:
                return True
    return False

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
    label_bam_file(bam, output + '.allelegroups.bam', file_to_read_to_allele_group['t'])
    if norm_bam is not None:
        label_bam_file(norm_bam, output + '.normal.allelegroups.bam', file_to_read_to_allele_group['n'])

def generate_new_vcf_header(in_vcf, norm_bam):
    with vcfpy.Reader.from_path(in_vcf) as reader:
        new_header = vcfpy.Header()
        for line in reader.header.lines:
            if type(line) in (vcfpy.HeaderLine, vcfpy.ContigHeaderLine):
                new_header.add_line(line)
    new_header.add_line(vcfpy.HeaderLine('source', 'FLAIR allelotyping'))
    new_header.add_line(vcfpy.InfoHeaderLine.from_mapping({'ID': 'AG', 'Number': 'R', 'Type': 'String', 'Description': 'Allele group(s) - can be considered a combo of phase set + haplotype'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'GT', 'Number': 1, 'Type': 'String', 'Description': 'Genotype'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'DP', 'Number': 1, 'Type': 'Integer', 'Description': 'Read depth'}))
    new_header.add_line(vcfpy.FormatHeaderLine.from_mapping({'ID': 'AD', 'Number': 'R', 'Type': 'Integer', 'Description': 'Allelic depths for the ref and alt alleles in the order listed'}))
    if norm_bam is not None:
        new_header.samples = vcfpy.SamplesInfos(['tumor', 'normal'])
    else:
        new_header.samples = vcfpy.SamplesInfos(['tumor'])
    return new_header

def write_vcf_file(output, new_header, variant_to_allele_group_counts_info, norm_bam, gene_to_allele_groups):
    allele_group_to_final_vars = {}
    for g in gene_to_allele_groups:
        for ag in gene_to_allele_groups[g]:
            allele_group_to_final_vars[ag] = set()
    with vcfpy.Writer.from_path(output + '.allelegroups.vcf', new_header) as writer:
        format_strings = ['GT', 'DP', 'AD']
        for (chrom, pos, ref, alt), var_data in variant_to_allele_group_counts_info.items():
            tot_cov = var_data['t_cov'] + var_data['n_cov']
            tot_var = var_data['t_var'] + var_data['n_var']
            # not reporting if homozygous
            # not reporting if in none or all allele groups
            if (tot_var != 0 and tot_var != tot_cov) and (len(var_data['ag']) != 0 and len(var_data['ag']) != gene_to_allele_groups[var_data['gene']]):
                # FIXME allow other types of alts
                alt_type = 'SNV'
                alt_desc = vcfpy.Substitution(alt_type, alt)
                these_calls = [vcfpy.Call('tumor', {'GT': '0/1', 'DP': var_data['t_cov'], 'AD': [var_data['t_cov'] - var_data['t_var'], var_data['t_var']]})]
                if norm_bam is not None:
                    these_calls.append(vcfpy.Call('normal', {'GT': '0/1', 'DP': var_data['n_cov'], 'AD': [var_data['n_cov'] - var_data['n_var'], var_data['n_var']]}))
                new_record = vcfpy.Record(chrom, int(pos), [], ref, [alt_desc], 100, ['PASS'], OrderedDict([('AG', var_data['ag'])]), format_strings, these_calls)
                writer.write_record(new_record)
                for allele_group in var_data['ag']:
                    # if allele_group not in allele_group_to_final_vars:
                    #     allele_group_to_final_vars[allele_group] = set()
                    allele_group_to_final_vars[allele_group].add((chrom, int(pos), ref, alt))
    return allele_group_to_final_vars

def write_allele_group_counts_read_map(index_to_allele_group_info, output, generate_map, norm_bam):
    with open(output + '.allelegroups.counts.tsv', 'w') as fh:
        fh_map = None
        if generate_map:
            fh_map = open(output + '.allelegroups.read.map.tsv', 'w')
        outline = ['#group_label', 'gene', 'tumor']
        if norm_bam is not None:
            outline.append('normal')
        fh.write('\t'.join(outline) + '\n')
        if generate_map:
            fh_map.write('\t'.join(outline) + '\n')
        for group_label, group_info in index_to_allele_group_info.items():
            ol = [group_label, group_info['gene'], group_info['t']]
            if norm_bam is not None:
                ol.append(group_info['n'])
            fh.write('\t'.join(ol[:2] + [str(len(x)) for x in ol[2:]]) + '\n')
            if generate_map:
                fh_map.write('\t'.join(ol[:2] + [','.join(x) for x in ol[2:]]) + '\n')
        if generate_map:
            fh_map.close()

def make_allele_group_label(group, allele_group_count, allele_group_to_reads, norm_bam, read_support):
    tot_reads_for_file = {'t': 0, 'n': 0}
    for file_label, read in allele_group_to_reads[group]:
        tot_reads_for_file[file_label] += 1
    allele_group_label = str(allele_group_count)
    # identify if is somatic
    if norm_bam is not None and tot_reads_for_file['n'] < read_support:
        allele_group_label += '-S'
    allele_group_count += 1
    return allele_group_count, allele_group_label


def get_allele_group_info(allele_group_to_reads, allele_group_count, var_info, gene_id, gene_chrom, index_to_allele_group_info, file_to_read_to_allele_group, variant_to_allele_group_counts_info, gene_to_allele_groups, norm_bam, read_support):
    gene_to_allele_groups[gene_id] = set()
    for group in allele_group_to_reads:
        allele_group_count, allele_group_label = make_allele_group_label(group, allele_group_count, allele_group_to_reads, norm_bam, read_support)
        gene_to_allele_groups[gene_id].add(allele_group_label)

        index_to_allele_group_info[allele_group_label] = {'gene': gene_id, 'vars': var_info, 'has_vars': group, 't': set(), 'n': set()}
        for file_label, read in allele_group_to_reads[group]:
            file_to_read_to_allele_group[file_label][read] = allele_group_label
            index_to_allele_group_info[allele_group_label][file_label].add(read)

        for varindex, is_var in group:
            pos, ref, alt = var_info[varindex].split(';')
            var_data = (gene_chrom, pos, ref, alt)
            if var_data not in variant_to_allele_group_counts_info:
                variant_to_allele_group_counts_info[var_data] = {'gene': gene_id, 'ag': [], 't_cov': 0, 't_var': 0, 'n_cov': 0, 'n_var': 0}
            if is_var == 1:
                variant_to_allele_group_counts_info[var_data]['ag'].append(allele_group_label)

    return allele_group_count

def get_variant_final_coverage(readvarinfo_to_reads, var_info, gene_chrom, variant_to_allele_group_counts_info):
    for og_group in readvarinfo_to_reads:
        for varindex, is_var in og_group:
            pos, ref, alt = var_info[varindex].split(';')
            var_data = (gene_chrom, pos, ref, alt)
            for file_label, read in readvarinfo_to_reads[og_group]:
                variant_to_allele_group_counts_info[var_data][file_label + '_cov'] += 1
                if is_var == 1:
                    variant_to_allele_group_counts_info[var_data][file_label + '_var'] += 1
