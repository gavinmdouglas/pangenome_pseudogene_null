#!/usr/bin/python3

from collections import defaultdict
from functions.io_utils import read_ids
import gzip
import os

# Get map of cd-hit-est input genes to specific clusters.
def main():

    species = read_ids('/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt')

    gene2sp = dict()
    gene2genome = dict()
    gene2contig = dict()

    for sp in species:
        sp_accessions = read_ids('/data1/gdouglas/projects/accessory_vs_pseudogene/genome_accessions/' + sp + '_accessions.txt')

        for acc in sp_accessions:

            pseudo_ids_filepath = '/data1/gdouglas/projects/accessory_vs_pseudogene/pseudogenes/processed/pseudo_filtered_ids/' + acc + '.txt'

            # Skip this accession if not in pseudogene output (meaning it was excluded).
            if not os.path.exists(pseudo_ids_filepath):
                continue

            pseudo_ids = read_ids(pseudo_ids_filepath, split_line = True)

            for pseudo_id in pseudo_ids:
                gene2genome[pseudo_id] = acc
                gene2sp[pseudo_id] = sp

            # Also read through the original pseudogene gff file to get the map of pseudogenes to contigs.
            pseudo_gff_file = '/data1/gdouglas/projects/accessory_vs_pseudogene/pseudogenes/pseudofinder_output/pseudos.gff/' + acc + '_pseudos.gff.gz'
            
            id_tally = 0

            with gzip.open(pseudo_gff_file, 'rt') as pseudo_gff_filehandle:
            
                for line in pseudo_gff_filehandle:
            
                    if line[0] == '#':
                        continue

                    line_split = line.split('\t')

                    contig = line_split[0]
                    pseudo_locus = line_split[8].split(';')[1].replace('locus_tag=', '').rstrip()

                    if pseudo_locus in pseudo_ids:
                        id_tally += 1

                        gene2contig[pseudo_locus] = contig

            if len(pseudo_ids) != id_tally:
                sys.exit('Error! Number of ids matched in gff does not match expected tally.')

            intact_gff_file = '/data1/gdouglas/projects/accessory_vs_pseudogene/pseudogenes/pseudofinder_output/intact.gff/' + acc + '_intact.gff.gz'
            
            with gzip.open(intact_gff_file, 'rt') as intact_gff_filehandle:
            
                for line in intact_gff_filehandle:
            
                    if line[0] == '#':
                        continue

                    line_split = line.split('\t')

                    contig = line_split[0]
                    intact_locus = line_split[8].replace('locus_tag=', '').rstrip()

                    gene2contig[intact_locus] = contig
                    gene2genome[intact_locus] = acc
                    gene2sp[intact_locus] = sp

    print("cluster\tgene\tcontig\taccession\tspecies")

    with open('/data1/gdouglas/projects/accessory_vs_pseudogene/cdhit_workflow_intact_and_intergenic/cdhit_filt_intergenic_pseudogenes_and_intact_c0.95_aL0.90_s0.90.clstr', 'r') as cdhitest_out:
        for line in cdhitest_out:
            line_split = line.split()

            if len(line_split) == 0:
                continue

            if line_split[0][0] == ">":
                current_cluster = 'c' + line_split[1]
                continue

            gene = line_split[2][1:]
            gene = gene.replace('.', '')
            print(current_cluster + '\t' + gene + '\t' + gene2contig[gene] + '\t' + gene2genome[gene] + '\t' + gene2sp[gene])


if __name__ == '__main__':
    main()
