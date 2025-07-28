#!/usr/bin/env python3
import sys
import argparse
from flair import FlairInputDataError
from flair.parse_iso_ids import IsoformInfo

def main():
    parser = argparse.ArgumentParser(description='options',
            usage='bed_to_gtf.py bed [options] > outfile.gtf')
    parser.add_argument('inputfile', type=str,
            action='store', help='isoforms in bed format')
    parser.add_argument('--force', action='store_true', dest='force',
            help='specify to not split isoform name by underscore into isoform and gene ids')
    args = parser.parse_args()

    bed_to_gtf(query=args.inputfile, force=args.force, outputfile='/dev/stdout')

def bed_to_gtf(query, outputfile, force=False):
    outfile = open(outputfile, 'w')
    for line in open(query):
        line = line.rstrip().split('\t')
        start = int(line[1])
        chrom, strand, score, name, start = line[0], line[5], line[4], line[3], int(line[1])
        tstarts = [int(n) + start for n in line[11].split(',')[:-1]]
        bsizes = [int(n) for n in line[10].split(',')[:-1]]
        end, thick_start, thick_end = int(line[2]), int(line[6]), int(line[7])

        if '|' not in name and not force:
            raise FlairInputDataError('Entry name should contain |-delimited transcriptid and geneid like so: \n'
                             'FL:4-1|ENST00000557334.5|ENSG00000133703.12 or FL:7-1|FL:7-1|chr12:25206000\n'
                             'So no GTF conversion was done. Please run identify_gene_isoform first\n'
                             'for best results, or run with --force')

        #if ';' in name:
        #    name = name.replace(';', ':')

        if force == True:
            transcript_id, gene_id = name, name
        else:
            isoinfo = IsoformInfo.parse_string(name)
            transcript_id, gene_id = isoinfo.full_transcript_name, isoinfo.gene_id

        attributes = 'gene_id \"{}\"; transcript_id \"{}\";'\
                                .format(gene_id, transcript_id)
        print('\t'.join([chrom, 'FLAIR', 'transcript', str(start+1), str(tstarts[-1]+bsizes[-1]), '.', strand, '.',
                attributes]), file=outfile)
        if thick_start != thick_end and (thick_start != start or thick_end != end):
            print('\t'.join([chrom, 'FLAIR', 'CDS', str(thick_start+1), str(thick_end), '.', strand, '.',
            attributes]), file=outfile)
            if strand == '+':
                print('\t'.join([chrom, 'FLAIR', 'start_codon', str(thick_start+1), str(thick_start+3), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '5UTR', str(start+1), str(thick_start+1), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '3UTR', str(thick_end), str(tstarts[-1]+bsizes[-1]), '.', strand, '.',
                attributes]), file=outfile)
            elif strand == '-':
                print('\t'.join([chrom, 'FLAIR', 'start_codon', str(thick_end-2), str(thick_end), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '3UTR', str(start+1), str(thick_start+1), '.', strand, '.',
                attributes]), file=outfile)
                print('\t'.join([chrom, 'FLAIR', '5UTR', str(thick_end), str(tstarts[-1]+bsizes[-1]), '.', strand, '.',
                attributes]), file=outfile)
        # if strand == '-':  # to list exons in 5'->3'
        #       for b in range(len(tstarts)):  # exon number
        #               bi = len(tstarts) - 1 - b  # block index
        #               attributes = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
        #                                               .format(gene_id, transcript_id, b)
        #               print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[bi]+1), \
        #                       str(tstarts[bi]+bsizes[bi]), '.', strand, '.', attributes]))
        # else:
        for b in range(len(tstarts)):
            attributes = 'gene_id \"{}\"; transcript_id \"{}\"; exon_number \"{}\";'\
                    .format(gene_id, transcript_id, b)
            print('\t'.join([chrom, 'FLAIR', 'exon', str(tstarts[b]+1),
                    str(tstarts[b]+bsizes[b]), '.', strand, '.', attributes]), file=outfile)

if __name__ == "__main__":
    main()
