#!/usr/bin/env python3
import sys
import csv
import os
from flair import FlairInputDataError
from flair.parse_iso_ids import IsoformInfo

def main():
    try:
        isoforms = sys.argv[1]
        s = float(sys.argv[2])
        outfilename = sys.argv[3]
    except:
        raise FlairInputDataError('usage: filter_isoforms_by_proportion_of_gene_expr.py isoforms support_percentage outfile')

    if s >= 1:
        raise FlairInputDataError('Support percentage should be a decimal e.g. 0.1 for 10%')

    filter_isoforms_by_proportion_of_gene_expr(isoforms=isoforms, outfilename=outfilename,
                                        support=s,)



def filter_isoforms_by_proportion_of_gene_expr(isoforms, outfilename, support):
    genes = dict()
    isoFH = open(isoforms, 'r')
    for line in isoFH:
        line = line.rstrip().split('\t')
        name = line[3]
        isoinfo = IsoformInfo.parse_string(name)
        gene = isoinfo.gene_id
        if gene not in genes:
            genes[gene] = set()
        genes[gene].add(tuple(line))

    with open(outfilename, 'wt') as outfile:
        writer = csv.writer(outfile, delimiter='\t', lineterminator=os.linesep)
        for gene in genes:
            gene_total = float(sum([float(iso[-1]) for iso in genes[gene]]))
            if gene_total == 0:
                continue
            for iso in genes[gene]:
                if float(iso[-1])/gene_total >= support:
                    writer.writerow(iso[:-1])

if __name__ == "__main__":
    main()
