#!/usr/bin/python3

import sys
import gzip
from functions.io_utils import read_ids
from collections import defaultdict

import rpy2.robjects as robjects
from rpy2.robjects import pandas2ri
pandas2ri.activate()
readRDS = robjects.r['readRDS']

# Get summary tables of which pseudogene and intact gene instances are redundant based on the presence of (another) intact gene of the same COG id in the same genome.
# This was originally written as an Rscript, but shifted to Python in hopes of speeding up the process.
# Requires input argument of species name to be given, also very much so hard-coded for use on Forsythia specifically.

def main():

    # Excluding soft.core as 'pseudogenes' in this partition are very likely false-positives and/or harder to interpret.
    partitions = ['ultra.cloud', 'other.cloud', 'shell']

    # Determine species being processed.
    all_species = read_ids('/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt')
    sp = sys.argv[1]
    if sp not in all_species:
        sys.exit('Stopping - input argument must be one of the ten input species names!')

    # Read in all COG-annotated clusters.
    COG_annot = pandas2ri.rpy2py(readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/all_filt_cluster_COG_annot.rds"))
    COG_annot_clusters = set(COG_annot.index)

    # Get sets of all species pseudogene/intact/both cluster ids and initialize breakdown objects.
    cluster_pseudogene_and_both_sets = dict()
    cluster_intact_and_both_sets = dict()

    pseudogene_breakdown = dict()
    intact_breakdown = dict()

    cluster_pangenome_categories_RDS = readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_pangenome_categories.rds")

    for partition in partitions:
        cluster_pseudogene_and_both_sets[partition] = set(cluster_pangenome_categories_RDS.rx2(sp).rx2("pseudo").rx2(partition)).union(set(cluster_pangenome_categories_RDS.rx2(sp).rx2("both").rx2(partition)))
        cluster_intact_and_both_sets[partition] = set(cluster_pangenome_categories_RDS.rx2(sp).rx2("intact").rx2(partition)).union(set(cluster_pangenome_categories_RDS.rx2(sp).rx2("both").rx2(partition)))

        pseudogene_breakdown[partition] = list()
        intact_breakdown[partition] = list()

    del cluster_pangenome_categories_RDS


    # Read in cluster ids that remained after length filters.
    length_filt_clusters = set()
    cluster_filt_tab = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz"
    with gzip.open(cluster_filt_tab, 'rt') as cluster_filt_tab_fh:
        cluster_filt_tab_fh.readline()
        for cluster_filt_line in cluster_filt_tab_fh:
            length_filt_clusters.add(cluster_filt_line.split()[0])

    # Note that the object for accessions with intact genes of a certain COG is a list.
    COG_intact_accessions_list = defaultdict(list)
    COG_pseudo_accessions = defaultdict(set)

    # Read through cluster membership breakdown table.
    cluster_member_breakdown_filename = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz"
    with gzip.open(cluster_member_breakdown_filename, "rt") as breakdown_fh:
        breakdown_header = breakdown_fh.readline().rstrip().split('\t')
        breakdown_colname_map = {index: element for element, index in enumerate(breakdown_header)}

        # Loop through all rows, and do the following:
        # 1. Keep only elements in this species and corresponding to COG-annotated clusters (and that passed filter check);
        # 2. Retain sets of accession ids that encode each intact COG (and also how many unique accessions do).
        # Also keep track of numbers of accessions that encode each pseudogene version of COG, so that these can be compared later;
        # 3. For all pseudogene lines, append them to separate lists of lists, depending on their partition type.

        for breakdown_line in breakdown_fh:
            breakdown_split = breakdown_line.rstrip().split('\t')

            cluster_id = breakdown_split[breakdown_colname_map['cluster']]

            if cluster_id not in length_filt_clusters:
                continue

            if breakdown_split[breakdown_colname_map['species']] != sp or cluster_id not in COG_annot_clusters:
                continue
            
            accession = breakdown_split[breakdown_colname_map['accession']]

            COG_genes = COG_annot.loc[cluster_id, 'COG']

            if '_pseudo' in breakdown_split[breakdown_colname_map['gene']]:

                for partition in partitions:
                    if cluster_id in cluster_pseudogene_and_both_sets[partition]:
                        pseudogene_breakdown[partition].append(tuple(breakdown_split[:-1] + [COG_genes]))
                        break

                # Get unique accession id that encodes each COG for later.
                for COG in COG_genes.split(','):
                    COG_pseudo_accessions[COG].add(accession)
            else:

                # Similar procedure for intact elements.
                for partition in partitions:
                    if cluster_id in cluster_intact_and_both_sets[partition]:
                        intact_breakdown[partition].append(tuple(breakdown_split[:-1] + [COG_genes]))
                        break

                for COG in COG_genes.split(','):
                    COG_intact_accessions_list[COG].append(accession)

    COG_intact_accessions = defaultdict(set)
    for COG in COG_intact_accessions_list.keys():
        COG_intact_accessions[COG] = set(COG_intact_accessions_list[COG])

    # First, write out table of numbers of accessions encode intact and pseudogenized versions of each COG gene family.
    tally_outfile = '/data1/gdouglas/projects/accessory_vs_pseudogene/summary/pseudo_compensation_breakdowns/' + sp + '_COG_accession_tallies.tsv'
    with open(tally_outfile, 'w') as tally_outfile_fh:
        print('COG\tintact\tpseudo', file = tally_outfile_fh)
        for COG in set(list(COG_intact_accessions.keys()) + list(COG_pseudo_accessions.keys())):
            print('\t'.join([COG,
                             str(len(COG_intact_accessions[COG])),
                             str(len(COG_pseudo_accessions[COG]))]),
                  file = tally_outfile_fh)


    # Then loop through all pseudogenes per partition and output potential compensation summary.
    redundant_outfile = '/data1/gdouglas/projects/accessory_vs_pseudogene/summary/pseudo_compensation_breakdowns/' + sp + '_redundant_COG_breakdown.tsv'
    with open(redundant_outfile, 'w') as redundant_outfile_fh:
        print('\t'.join(['element', 'cluster', 'partition', 'accession', 'COG',  'redundant_intact_COG']),
              file = redundant_outfile_fh)

        for partition in partitions:

            # For pseudogene elements.
            for info in pseudogene_breakdown[partition]:

                COG_genes = info[-1]

                accession = info[breakdown_colname_map['accession']]

                redundant = 'FALSE'
                for COG in COG_genes.split(','):
                    if accession in COG_intact_accessions[COG]:
                        redundant = 'TRUE'
                        break

                print('\t'.join([info[breakdown_colname_map['gene']],
                                 info[breakdown_colname_map['cluster']],
                                 partition,
                                 accession,
                                 COG_genes,
                                 redundant]),
                      file = redundant_outfile_fh)

            # For intact elements.
            for intact_info in intact_breakdown[partition]:

                COG_genes = intact_info[-1]

                accession = intact_info[breakdown_colname_map['accession']]

                redundant = 'FALSE'
                for COG in COG_genes.split(','):
                    if COG_intact_accessions_list[COG].count(accession) >= 2:
                        redundant = 'TRUE'
                        break

                print('\t'.join([intact_info[breakdown_colname_map['gene']],
                                 intact_info[breakdown_colname_map['cluster']],
                                 partition,
                                 accession,
                                 COG_genes,
                                 redundant]),
                      file = redundant_outfile_fh)


if __name__ == '__main__':
    main()
