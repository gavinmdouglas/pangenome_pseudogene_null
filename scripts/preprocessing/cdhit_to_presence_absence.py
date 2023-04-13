#!/usr/bin/python3

import argparse
from collections import defaultdict


def cluster_presence_info(genomes, cluster_info):
    out = []
    for g in genomes:
        if g in cluster_info.keys():
            out.append(','.join(sorted(cluster_info[g])))
        else:
            out.append('')
    return(out)


def main():

    parser = argparse.ArgumentParser(

        description='Given a cd-hit cluster output file, and a mapfile of sequence element ids to genome accessions, '
                    'return table giving breakdown of cluster presence/absence across accessions.',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--cdhit", metavar="PATH", type=str, required=True,
                        help="cd-hit cluster outputfile")

    parser.add_argument("--map", metavar="FILE", type=str, required=True,
                        help="Whitespace-delimited mapfile of element ids to genome accessions")

    parser.add_argument("--prefix", metavar="STRING", type=str, required=False, default='',
                        help="Prefix to add to cluster output names.")

    parser.add_argument("--rm_first_field", action='store_true',
                        help="Set to remove first field in sequence element ids after splitting by _ in cluster output file.",
                        required=False, default=False)

    args = parser.parse_args()

    # Read in element ids to accessions
    global element_to_accession
    element_to_accession = dict()
    accessions = set()
    with open(args.map, 'r') as map_fh:
        for map_line in map_fh:
            map_line_split = map_line.split()
            accessions.add(map_line_split[1])
            element_to_accession[map_line_split[0]] = map_line_split[1]

    # Read through cd-hit cluster output file, and output prevalence info across accessions.
    current = None
    cluster_set = defaultdict(list)

    accessions_sorted = sorted(list(accessions))
    print('\t'.join(['cluster'] + accessions_sorted))

    with open(args.cdhit, 'r') as cdhit_fh:
        for cdhit_line in cdhit_fh:
            if cdhit_line.startswith('>Cluster '):
                if current:
                    output_line_split = [current] + \
                                        cluster_presence_info(genomes=accessions_sorted,
                                                              cluster_info=cluster_set)

                    print('\t'.join(output_line_split))

                    current = None
                    cluster_set = defaultdict(list)

                current = args.prefix + cdhit_line[1:].rstrip().replace('Cluster ', 'c')
            else:
                cdhit_line_split = cdhit_line.split()
                element = cdhit_line_split[2][1:-3]

                if args.rm_first_field:
                    element_split = element.split('_')
                    element = '_'.join(element_split[1:])

                element_accession = element_to_accession[element]
                cluster_set[element_accession].append(element)

    if current:
        output_line_split = [current] + \
                            cluster_presence_info(genomes=accessions_sorted,
                                                  cluster_info=cluster_set)

        print('\t'.join(output_line_split))


if __name__ == '__main__':
    main()
