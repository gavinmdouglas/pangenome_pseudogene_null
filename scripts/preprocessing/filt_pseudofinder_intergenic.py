#!/usr/bin/python3

import sys
import argparse
import gzip
from collections import defaultdict


def main():

    parser = argparse.ArgumentParser(

        description='Get filtered intergenic pseudogene ids per species. '
                    'Read through pseudofinder-output pseudogene GFFs given folder name containing all GFFs and BLASTX files and file containing genome ids. '
                    'First, restrict to intergenic pseudogenes only. And make sure that these pseudogenes are at least 500 bp away from contig edge. '
                    'Then, only keep intergenic pseudogenes that are between (inclusively) 100-5000 bp in length. '
                    'Last, compare pseudogene length with mean size of database sequence hits (mean ref. - pseudogene). '
                    'Only keep pseudogenes with size differences inclusively between -500 and 2000.',

        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--pseudofinder_outfolder", metavar="FOLDER", type=str, required=True,
                        help="Location of all pseudofinder outfiles (including GFFs and BLASTX hits)")

    parser.add_argument("--species_ids", metavar="FILE", type=str, required=True,
                        help="File with genome accessions for species: one per line (assumed to be prefix in pseudofinder outfiles).")

    parser.add_argument("--database_lengths", metavar="FILE", type=str, required=True,
                        help="Gzipped, tab-delimited file (with header) with nucleotide length per database sequence.")

    parser.add_argument("--stats_outfile", metavar="OUT", type=str, required=True,
                        help="File to write filtering stats per genome.")

    args = parser.parse_args()

    # Read in genome ids for species.
    genome_ids = set()
    with open(args.species_ids, 'r') as species_ids_fh:
        for species_id_line in species_ids_fh:
            genome_ids.add(species_id_line.rstrip())

    # First read through BLASTX files once to get all database
    # sequences for which the length needs to be parsed.
    database_hits = set()
    for genome_id in genome_ids:
        blastx_file = args.pseudofinder_outfolder + '/' + genome_id + '_intergenic.fasta.blastX_output.tsv'
        with open(blastx_file, 'r') as blastx_file_fh1:
            for blastx_line1 in blastx_file_fh1:
                database_hits.add(blastx_line1.split()[1])

    # Then read through actual database length file and get
    # nucleotide lengths for the subset of hits identified.
    ref_lengths = dict()
    with gzip.open(args.database_lengths, mode="rt") as db_lengths_fh:
        db_lengths_fh.readline()
        for db_line in db_lengths_fh:
            db_line_split = db_line.split()
            if db_line_split[0] in database_hits:
                ref_lengths[db_line_split[0]] = int(db_line_split[1])

    # Then read through the BLASTX files again and get the mean size of all reference hits.
    summed_sizes = defaultdict(int)
    num_hits = defaultdict(int)
    for genome_id in genome_ids:
        blastx_file = args.pseudofinder_outfolder + '/' + genome_id + '_intergenic.fasta.blastX_output.tsv'
        with open(blastx_file, 'r') as blastx_file_fh2:
            for blastx_line2 in blastx_file_fh2:
                blastx_line_split = blastx_line2.split()
                ign_id = blastx_line_split[0]
                ref_id = blastx_line_split[1]
                summed_sizes[ign_id] += ref_lengths[ref_id]
                num_hits[ign_id] += 1

    mean_sizes = dict()
    for intergenic_id in summed_sizes.keys():
        mean_sizes[intergenic_id] = summed_sizes[intergenic_id] / num_hits[intergenic_id]

    # Then parse all gff files and determine set of intergenic pseudogenes to retain.
    # Keep track of total count failed due to different reasons too.
    stats_fh = open(args.stats_outfile, 'w')
    print('\t'.join(['genome', 'passed', 'raw', 'edge_fail',
                     'size_fail', 'ref_diff_fail']),
          file=stats_fh)

    for genome_id in genome_ids:
        raw = 0
        edge_fail = 0
        size_fail = 0
        ref_diff_fail = 0
        passed = 0
        gff_file = args.pseudofinder_outfolder + '/' + genome_id + '_pseudos.gff'
        contig_sizes = dict()
        with open(gff_file, 'r') as gff_file_fh:
            for gff_line in gff_file_fh:

                # Parse contig sizes.
                if gff_line.startswith('##sequence-region'):
                    gff_split = gff_line.split()
                    contig_sizes[gff_split[1]] = int(gff_split[3])
                    continue

                # Only consider intergenic pseudogenes.
                if 'Reason(s):Intergenic region ' not in gff_line:
                    continue

                # Keep tally of raw intergenic pseudogenes parsed.
                raw += 1

                line_split = gff_line.split('\t')

                # Get pseudogene id.
                if ';old_locus_tag=' not in line_split[8]:
                    sys.exit('Stopping - no old locus tag!')
                pseudo_id_split = line_split[8].split(';')
                old_locus_tag = pseudo_id_split[-1].replace('old_locus_tag=', '').rstrip()
                pseudo_id = pseudo_id_split[-2].replace('locus_tag=', '')

                dist_from_contig_end = contig_sizes[line_split[0]] - int(line_split[4])
                pseudogene_start = int(line_split[3])
                failed = False

                # Check that hit is not too close to contig end.
                if dist_from_contig_end < 500 or pseudogene_start < 500:
                    failed = True
                    edge_fail += 1

                # Check that pseudogene is not too small or too large.
                pseudogene_size = int(line_split[4]) - int(line_split[3]) + 1
                if pseudogene_size < 100 or pseudogene_size > 5000:
                    failed = True
                    size_fail += 1

                # Check that pseudogene is not too different in size from mean reference hit.
                diff_vs_ref = mean_sizes[old_locus_tag] - pseudogene_size
                if diff_vs_ref < -500 or diff_vs_ref > 2000:
                    failed = True
                    ref_diff_fail += 1

                # Then print out pseudogene id if passed.
                if not failed:
                    passed += 1
                    print(pseudo_id + '\t' + genome_id)

        # Write out stats per genome.
        print('\t'.join([genome_id, str(passed), str(raw), str(edge_fail),
                         str(size_fail), str(ref_diff_fail)]),
              file=stats_fh)

    stats_fh.close()


if __name__ == '__main__':
    main()
