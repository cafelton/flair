#! /usr/bin/env python3

import os
import sys
import argparse
import pipettor
import pysam
from shutil import rmtree
import logging
from flair.io_utils import make_temp_dir
from flair.pycbio.hgdata.bed import BedReader
from flair.flair_bed import FlairBed
from flair.flair_transcriptome import _check_junction_subset
from flair.read_processing import generate_genomic_alignment_read_to_clipping_file
from flair.count_sam_transcripts import run_count_sam_transcripts
from flair.counts_to_tpm import counts_to_tpm
import multiprocessing as mp


os.environ['OPENBLAS_NUM_THREADS'] = '1'

def parse_args():
    parser = argparse.ArgumentParser()
    required = parser.add_argument_group('required named arguments')
    required.add_argument('--manifest', action='store', type=str,
                          required=True, help='Tab delimited file containing sample id, condition, batch, reads.fq')
    required.add_argument('--isoform_fa', action='store',
                          type=str, required=True, help='FastA of FLAIR collapsed isoforms')
    parser.add_argument('-o', '--output', type=str, action='store', default='flair.quantify',
                        help='''output file name base for FLAIR quantify (default: flair.quantify)''')
    parser.add_argument('-t', '--threads', type=int,
                        action='store', default=4, help='minimap2 number of threads (4)')
    parser.add_argument('--sample_id_only', default=False, action='store_true',
                        help='''only use sample id in output header''')
    parser.add_argument('--tpm', action='store_true', default=False,
                        help='Convert counts matrix to transcripts per million and output as a separate file named <output>.tpm.tsv')
    parser.add_argument('--quality', type=int, action='store', default=0,
                        help='''minimum MAPQ of read assignment to an isoform (0)''')
    parser.add_argument('--trust_ends', default=False, action='store_true',
                        help='specify if reads are generated from a long read method with minimal fragmentation')
    parser.add_argument('--generate_map', default=False, action='store_true',
                        help='''create read-to-isoform assignment files for each sample (default: not specified)''')
    required.add_argument('--isoform_bed', required=True, type=str, action='store',
                          help='''isoform .bed file''')
    # parser.add_argument('--stringent', default=False, action='store_true',
    #                     help='''Supporting reads must cover 80 percent of their isoform and extend at least 25 nt into the
    #                     first and last exons. If those exons are themselves shorter than 25 nt, the requirement becomes
    #                     'must start within 4 nt from the start" or "must end within 4 nt from the end" ''')
    # parser.add_argument('--check_splice', default=False, action='store_true',
    #                     help='''enforce coverage of 4 out of 6 bp around each splice site and no
    #                     insertions greater than 3 bp at the splice site''')
    # parser.add_argument('--output_bam', default=False, action='store_true',
    #                     help='[for development] whether to output bam file of reads aligned to correct isoforms')
    args = parser.parse_args()
    return args

def check_args(args):
    args = parse_args()
    # if (args.stringent or args.check_splice):
    #     if not args.isoforms:
    #         raise Exception('Please specify isoform models as .bed file using --isoform_bed')
    if not os.path.exists(args.isoform_bed):
        raise Exception('Isoform models bed file path does not exist: ' + args.isoform_bed)
    elif args.isoform_bed.endswith('.psl'):
        raise Exception('** Error. Flair no longer accepts PSL input. Please use psl_to_bed first.')
    if not os.path.exists(args.isoform_fa):
        raise Exception('Isoform sequences fasta file path does not exist: ' + args.isoform_fa)
    if not os.path.exists(args.manifest):
        raise Exception('Manifest file path does not exist: ' + args.manifest)

def load_manifest(manifest, sample_id_only):
    sample_data = []
    for line in open(manifest):
        cols = line.rstrip().split('\t')
        if len(cols) < 4:
            raise Exception(f'Expected 4 columns in tab-delimited manifest.tsv, got {len(cols)}. Exiting.')

        sample, group, batch, readFile = cols
        if sample_id_only is False:
            if '_' in sample or '_' in group or '_' in batch:
                raise Exception(f'Please do not use underscores in the id, condition, or batch fields of {args.r}.')

        if not os.path.exists(readFile):
            raise Exception('Query file path does not exist: {}'.format(readFile))

        sample_data.append((sample, group, batch, readFile))
    return sample_data

class GeneData():
    def __init__(self, gene_id, chrom, strand):
        self.gene_id = gene_id
        self.chrom = chrom
        self.strand = strand
        self.left_bound = None
        self.right_bound = None
        self.isoform_beds = []
        self.isoform_fas = []

    def add_isoform_bed(self, isoform):
        self.isoform_beds.append(isoform)
        if self.left_bound is None or isoform.chromStart < self.left_bound:
            self.left_bound = isoform.chromStart
        if self.right_bound is None or isoform.chromEnd > self.right_bound:
            self.right_bound = isoform.chromEnd

def load_isoform_data(isoform_bed, isoform_fa):
    gene_data = {}
    isoform_to_gene = {}
    for bed in BedReader(isoform_bed, bedClass=FlairBed):
        if bed.gene_id not in gene_data:
            gene_data[bed.gene_id] = GeneData(bed.gene_id, bed.chrom, bed.strand)
        gene_data[bed.gene_id].add_isoform_bed(bed)
        isoform_to_gene[bed.name] = bed.gene_id
    iso_name = None
    for line in open(isoform_fa):
        if line[0] == '>':
            iso_name = line[1:].rstrip().split()[0]
        else:
            gene_id = isoform_to_gene[iso_name]
            gene_data[gene_id].isoform_fas.append((iso_name, line.rstrip()))
    return gene_data

def write_unique_bound(fh, isoform, unique_seq_bound):
    unique_seq_bound = list(set(unique_seq_bound))
    if isoform.strand == '-':
        # just invert the indexes
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{abs(unique_seq_bound[i][0] - 1)}_{unique_seq_bound[i][1]}'
    else:
        for i in range(len(unique_seq_bound)):
            unique_seq_bound[i] = f'{unique_seq_bound[i][0]}_{unique_seq_bound[i][1]}'
    if len(unique_seq_bound) > 0:
        fh.write(isoform.name + '\t' + ','.join(unique_seq_bound) + '\n')

def load_unique_bound(temp_prefix, gene_info):
    with open(temp_prefix + 'isoforms.uniquebound.txt', 'w') as fh:
        for isoform in gene_info.isoform_beds:
            if len(isoform.blocks) > 1:
                unique_seq_bound = []
                terminal_exon_is_subset = [0, 0]  # first exon is a subset, last exon is a subset
                superset_support = []
                juncs = isoform.getGaps()
                first_exon, last_exon = isoform.blocks[0], isoform.blocks[-1]
                for otheriso in gene_info.isoform_beds:
                    if isoform.name != otheriso.name and len(otheriso.blocks) > 1:
                        _check_junction_subset(juncs, first_exon, last_exon, otheriso.read_support, otheriso.getGaps(), otheriso.blocks,
                                               terminal_exon_is_subset, superset_support, unique_seq_bound)
                write_unique_bound(fh, isoform, unique_seq_bound)             

def write_combined_counts(sample_data, gene_data, temp_dir, output, sample_id_only):
    with open(output + '.counts.tsv', 'w') as fh:
        if sample_id_only:
            fh.write('\t'.join(['ids'] + [x[0] for x in sample_data]) + '\n')
        else:
            fh.write('\t'.join(['ids'] + ['_'.join(x[:3]) for x in sample_data]) + '\n')
        for gene_id in gene_data:
            gene_info = gene_data[gene_id]
            iso_to_counts = {x.name: [0] * len(sample_data) for x in gene_info.isoform_beds}
            for i, (sample, group, batch, bamfile) in enumerate(sample_data):
                for line in open(temp_dir + gene_id + '/' + sample + '.isoform.counts.txt'):
                    line = line.rstrip().split('\t')
                    iso_to_counts[line[0]][i] = int(line[1])
            for iso in iso_to_counts:
                fh.write('\t'.join([iso + '_' + gene_id] + [str(x) for x in iso_to_counts[iso]]) + '\n')

def write_map_out(sample_data, gene_data, temp_dir, output, generate_map):
    if generate_map:
        for sample, group, batch, bamfile in sample_data:
            with open(f'{output}.{sample}.read.map.txt', 'w') as fh:
                for gene_id in gene_data:
                    for line in open(temp_dir + gene_id + '/' + sample + '.isoform.read.map.txt'):
                        fh.write(line)

def get_counts_for_sample(sample, bamfile, temp_prefix, gene_info, generate_map, trust_ends):
    temp_prefix_sample = temp_prefix + sample
    pipettor.run([('samtools', 'view', '-h', bamfile, gene_info.chrom + ':' + str(gene_info.left_bound) + '-' + str(gene_info.right_bound)),
                  ('samtools', 'fasta', '-')],
                 stdout=temp_prefix_sample + '.reads.fasta')
    bam_file = pysam.AlignmentFile(bamfile, 'rb')
    generate_genomic_alignment_read_to_clipping_file(temp_prefix_sample, bam_file, gene_info.chrom, gene_info.left_bound, gene_info.right_bound)
    read_map_file = temp_prefix_sample + '.isoform.read.map.txt' if generate_map else None
    mm2_cmd = ['minimap2', '-a', '-N', '4', '--MD', temp_prefix + 'isoforms.fa', temp_prefix_sample + '.reads.fasta']
    run_count_sam_transcripts(
        mm2_cmd=mm2_cmd,
        output=temp_prefix_sample + '.isoform.counts.txt',
        trimmedreads=temp_prefix_sample + '.reads.genomicclipping.txt',
        generate_map=read_map_file,
        end_norm_dist=0,
        stringent=True,
        allow_UTR_indels=True,  # is_annot,
        output_bam=False,  # args.output_bam,
        check_splice=True,
        isoforms=temp_prefix + 'isoforms.bed',
        trust_ends=trust_ends,
        unique_bound=temp_prefix + 'isoforms.uniquebound.txt',
    )

def get_counts_for_gene(input):
    temp_dir, gene_id, gene_info, sample_data, generate_map, trust_ends = input
    logging.debug(f'realigning reads to {gene_id}')
    temp_prefix = temp_dir + gene_id + '/'
    os.makedirs(temp_prefix)
    with open(temp_prefix + 'isoforms.bed', 'w') as fh:
        for bed in gene_info.isoform_beds:
            bed.write(fh)
    with open(temp_prefix + 'isoforms.fa', 'w') as fh:
        for name, seq in gene_info.isoform_fas:
            fh.write('>' + name + '\n' + seq + '\n')
    load_unique_bound(temp_prefix, gene_info)

    for sample, group, batch, bamfile in sample_data:
        get_counts_for_sample(sample, bamfile, temp_prefix, gene_info, generate_map, trust_ends)

def quantify():
    logging.info('loading isos')
    args = parse_args()
    check_args(args)
    temp_dir = make_temp_dir(args.output)
    sample_data = load_manifest(args.manifest, args.sample_id_only)
    gene_data = load_isoform_data(args.isoform_bed, args.isoform_fa)

    logging.info(f'Re-aligning reads to transcriptome. Writing temporary files into {temp_dir}')

    packed = []
    for gene_id in gene_data:
        packed.append((temp_dir, gene_id, gene_data[gene_id], sample_data, args.generate_map, args.trust_ends))

    if args.threads == 1:
        for p in packed:
            get_counts_for_gene(p)
    else:
        mp.set_start_method('fork', force=True)
        with mp.Pool(args.threads) as pool:
            pool.map(get_counts_for_gene, packed)

    logging.info('writing quantify output')
    write_combined_counts(sample_data, gene_data, temp_dir, args.output, args.sample_id_only)

    write_map_out(sample_data, gene_data, temp_dir, args.output, args.generate_map)

    counts_to_tpm(open(args.output + '.counts.tsv'), args.output + '.tpm.tsv')

    rmtree(temp_dir)


if __name__ == '__main__':
    # FIXME: need proper error handling
    sys.exit(quantify())
