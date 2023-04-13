#!/usr/bin/python3

import argparse
import sys
import os
import textwrap

from functions.io_utils import read_fasta, read_ids


def main():

    parser = argparse.ArgumentParser(

    description="Reads in FASTA and will only write out sequences "
                "with ids in the specified file. Note that the sequence "
                "id is taken to be the first field after splitting by "
                "whitespace and removing the > character.",

    formatter_class=argparse.RawDescriptionHelpFormatter)


    parser.add_argument("-f", "--fasta", metavar="FASTA", type=str,
                        required=True, help="Path to FASTA file")

    parser.add_argument("--ids2keep", metavar="FILE", type=str, required=True,
                        help="File with ids to keep - one per line. If there are multiple whitespace delimited columns, it will take the first column as the id.")

    parser.add_argument("-o", "--output", metavar="OUTPUT", type=str,
                        required=True, help="Path to output FASTA file")

    args = parser.parse_args()

    seqs = read_fasta(args.fasta, cut_header=True)
    
    ids2keep = read_ids(filename = args.ids2keep,
                        split_line = True)

    num_kept = 0

    out_fasta = open(args.output, 'w')

    for id2keep in ids2keep:
        if id2keep not in seqs.keys():
            sys.exit('Stopping, id ' + id2keep + ' not found in input FASTA.')

        num_kept += 1

        out_fasta.write(">" + id2keep + "\n")
        out_fasta.write(textwrap.fill(seqs[id2keep], width=70) + "\n")

    out_fasta.close()


    if num_kept != len(ids2keep):
        sys.exit('Error - number of seqs printed out does not match expected ids to keep specified. Kept: ' + str(num_kept))


if __name__ == '__main__':
    main()

