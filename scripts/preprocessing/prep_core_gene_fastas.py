#!/usr/bin/python3

import sys
import argparse
from functions.io_utils import write_fasta


# Get FASTA files of raw sequences for all core genes in panaroo output.
def main():

    parser = argparse.ArgumentParser(

    description="Reads in panaroo gene_data.csv and gene_presence_absence.csv files and output raw FASTA per core gene. Not that only genes without the flags of refound, stop, and/or len are retained.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


    parser.add_argument("-d", "--data", metavar="Gene-data", type=str,
                        required=True, help="Path to gene_data.csv.")

    parser.add_argument("-p", "--presence", metavar="Gene-data", type=str,
                        required=True, help="Path to gene_presence_absence.csv.")

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        required=True, help="Path to folder to output FASTAs.")

    args = parser.parse_args()

    gene_seqs = dict()
    col_map = dict()
    header_flag = True

    with open(args.data, 'r') as data_fh:
        for data_line in data_fh:
            data_split = data_line.split(',')

            if header_flag:
                for i in range(len(data_split)):
                    col_map[data_split[i]] = i
                header_flag = False
                continue

            annot_i = col_map['annotation_id']
            dna_i = col_map['dna_sequence']
            gene_seqs[data_split[annot_i]] = data_split[dna_i]


    with open(args.presence, 'r') as presence_fh:
        next(presence_fh)
        for presence_line in presence_fh:
            presence_line = presence_line.rstrip()
            presence_split = presence_line.split(',')
            name = presence_split[0]
            genes = presence_split[3:]

            passed = True
            tmp_seqs = dict()
            for g in genes:
                if g == '' or ';' in g or '_len' in g or '_stop' in g or '_refound' in g:
                    passed = False
                    break

                tmp_seqs[g] = gene_seqs[g]

            if passed:
                outfile = args.output + '/' + name + '.fna'
                write_fasta(seq = tmp_seqs, outfile = outfile)


if __name__ == '__main__':

    main()
