#! /usr/bin/env python3

import argparse
import pysam
from gtf_io import load_gtf_to_gene_data, GtfExon
from copy import deepcopy
import graphviz
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
import flair.variant_processing as vp


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
        parser.error(f'At least one of gtf and isoform_bed must be provided')
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

def check_read_is_in_gene(gene_juncs, introns):
    if len(set(gene_juncs) & set(introns)) > 0:
        return True
    return False

def load_bam_file_for_region(bam, bam_label, rchrom, rstart, rend, my_vars, genome, gene_juncs, readvarinfo_to_reads):
    with pysam.AlignmentFile(bam, 'rb') as bam:
        for a in bam.fetch(rchrom, rstart, rend):
            if a.is_mapped and not a.is_secondary and not a.is_supplementary:
                introns, insertions, deletions = vp.parse_cigar_for_introns_indels(a.cigartuples, a.reference_start)
                if check_read_is_in_gene(gene_juncs, introns):
                    read_id, read_vars = vp.load_bam_read_vars(a, bam_label, deepcopy(my_vars), genome, insertions, deletions)
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
    vartoalt = vp.combine_vcf_files(my_vcfs)
    chrom_to_region_to_vars = vp.reorganize_vars(vartoalt)
    gene_to_all_exons, gene_to_all_juncs = load_isoform_data(args.gtf, args.isoform_bed)

    gene_to_vars = vp.get_gene_to_all_vars(gene_to_all_exons, chrom_to_region_to_vars)
    # tempdir = make_temp_dir(args.output)

    print('loading bam files and creating allele groups')
    file_to_read_to_allele_group = {'t': {}, 'n': {}}
    index_to_allele_group_info = {}
    variant_to_allele_group_counts_info = {}
    gene_to_allele_groups = {}
    allele_group_count = 1
    # TODO can parallelize at the gene level, will need to write to intermediate files
    for gene_id in gene_to_vars:
        # only process reads if gene has variants
        # and gene is spliced (junctions are how we assign reads to gene)
        if len(gene_to_vars[gene_id]) > 0 and len(gene_to_all_juncs[gene_id]) > 0:
            gene_chrom, gene_start, gene_end = vp.get_gene_boundaries(gene_to_all_exons[gene_id])
            my_vars = vp.create_gene_vars_dict(gene_to_vars[gene_id])
            var_info, _ = vp.simplify_gene_vars(my_vars)
            readvarinfo_to_reads = load_bam_files_for_region(args.bam, args.norm_bam, gene_chrom, gene_start, gene_end, my_vars, genome, gene_to_all_juncs[gene_id])
            allele_group_to_reads = vp.process_alleotype_graph(deepcopy(readvarinfo_to_reads), args.read_support, args.frac_support)

            # don't report if there's only one final group
            if len(allele_group_to_reads) > 1:
                allele_group_count = vp.get_allele_group_info(allele_group_to_reads, allele_group_count, var_info, gene_id, gene_chrom, index_to_allele_group_info,
                                                              file_to_read_to_allele_group, variant_to_allele_group_counts_info, gene_to_allele_groups, args.norm_bam, args.read_support)

                # getting variant coverage from original assignments, not collapsed ones
                vp.get_variant_final_coverage(readvarinfo_to_reads, var_info, gene_chrom, variant_to_allele_group_counts_info)

    print('writing vcf file')
    new_header = vp.generate_new_vcf_header(args.vcf, args.norm_bam)
    vp.write_vcf_file(args.output, new_header, variant_to_allele_group_counts_info, args.norm_bam, gene_to_allele_groups)

    vp.write_allele_group_counts_read_map(index_to_allele_group_info, args.output, args.generate_map, args.norm_bam)

    if args.annotate_bams:
        print('labeling bam files')
        vp.label_bam_files(args.bam, args.norm_bam, args.output, file_to_read_to_allele_group)


if __name__ == "__main__":
    getvariants()
