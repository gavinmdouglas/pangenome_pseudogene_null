#!/usr/bin/python3

import sys
import gzip
from functions.io_utils import read_ids
from collections import defaultdict
import numpy

# Get tab-delimited table of intergenic pseudogene ids (filtered based on kept accessions + within size cutoff) to UniRef matches.

def main():

    old2new_ids = dict()
    # Read in filtered pseudogene ids.
    with gzip.open("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_pseudogene_ids_acc_and_length_filt.tsv.gz", mode="rt") as filt_pseudo_ids:
        filt_pseudo_ids.readline()
        for line in filt_pseudo_ids:
            line = line.rstrip()
            line_split = line.split('\t')
            old2new_ids[line_split[1]] = line_split[0]

    species = read_ids('/data1/gdouglas/projects/hgt_fragments/mapfiles/species.txt')

    accessions = set()
    for sp in species:
        sp_acc = read_ids('/data1/gdouglas/projects/hgt_fragments/genome_accessions/' + sp + '_accessions.txt')
        for acc in sp_acc:
            accessions.add(acc)

    print('intergenic_pseudo_id\traw_id\tuniref_hit')

    for acc in accessions:

        intergenic_pseudo_blastx = '/data1/gdouglas/projects/hgt_fragments/pseudogenes/pseudofinder_output/intergenic.fasta.blastX_output/' + acc + '_intergenic.fasta.blastX_output.tsv.gz'

        # Read through BLASTX-formatted output.
        # These are the columns:
        # qseqid sseqid pident slen mismatch gapopen qstart qend sstart send evalue bitscore stitle
        with gzip.open(intergenic_pseudo_blastx, mode="rt") as intergenic_pseudo_blastx_in:

            for intergenic_blastx_line in intergenic_pseudo_blastx_in:

                line_split = intergenic_blastx_line.split('\t')

                if line_split[0] not in old2new_ids.keys():
                    continue

                print('\t'.join([old2new_ids[line_split[0]], line_split[0], line_split[1]]))


if __name__ == '__main__':
    main()
