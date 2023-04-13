#!/usr/bin/python3

import argparse
import sys
import os
import textwrap
from functions.io_utils import read_fasta


def main():

    parser = argparse.ArgumentParser(
        description="Subset pseudogene FASTAs based on file of ids, for species separately (based on genome ids file).",
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--pseudofinder_outfolder", metavar="FOLDER", type=str, required=True,
                        help="Location of all pseudofinder outfiles (including GFFs and BLASTX hits)")

    parser.add_argument("--species_ids", metavar="FILE", type=str, required=True,
                        help="File with genome accessions for species: one per line (assumed to be prefix in pseudofinder outfiles).")

    parser.add_argument("--ids2keep", metavar="FILE", type=str, required=True,
                        help="File with pseudogene ids to keep - one per line (first field).")

    parser.add_argument("--sp_full2short", metavar="FILE", type=str, required=True,
                        help="Mapfile of species names to short ids (the short form will be prefixed to pseudogene ids)")

    args = parser.parse_args()

    genome_ids = set()
    with open(args.species_ids, 'r') as genome_file:
        for genome_line in genome_file:
            genome_line_split = genome_line.split()
            if len(genome_line_split) > 0:
                genome_ids.add(genome_line_split[0])

    ids2keep = set()
    with open(args.ids2keep, 'r') as ids_file:
        for line in ids_file:
            line_split = line.split()
            if len(line_split) > 0:
                ids2keep.add(line_split[0])

    # Assume that species name is the species id filename
    # after removing the extension.
    species, file_extension = os.path.splitext(args.species_ids)
    species = os.path.basename(species)

    # Figure out the short-form for species name.
    with open(args.sp_full2short, 'r') as map_file:
        for map_line in map_file:
            map_file_split = map_line.split()
            if map_file_split[0] == species:
                species_short = map_file_split[1]
                break

    # Loop through all species genome ids and read in FASTAs of pseudogenes.
    # This is assuming that this will only use a little memory.
    all_pseudogenes = dict()

    for genome_id in genome_ids:
        fasta_file = args.pseudofinder_outfolder + '/' + genome_id + '_pseudos.fasta'
        genome_pseudo = read_fasta(fasta_file, cut_header=True)

        if len(all_pseudogenes.keys() & genome_pseudo.keys()) > 0:
            sys.exit('Error - duplicate pseudogene ids across genomes.')

        all_pseudogenes.update(genome_pseudo)

    # Then print out pseudogenes.
    for id2keep in ids2keep:
        if id2keep not in all_pseudogenes.keys():
            sys.exit('Stopping, id ' + id2keep + ' not found in input FASTAs.')

        print('>' + species_short + '_' + id2keep)
        print(textwrap.fill(all_pseudogenes[id2keep], width=70))


if __name__ == '__main__':
    main()
