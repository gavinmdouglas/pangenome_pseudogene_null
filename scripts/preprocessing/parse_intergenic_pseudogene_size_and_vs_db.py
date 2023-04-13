#!/usr/bin/python3

import sys
import gzip
from functions.io_utils import read_ids
from collections import defaultdict
import numpy

# Get absolute size of intergenic pseudogenes, and also average/SD percent size versus database matches.

def main():

    uniref90_lengths = dict()
    # Read in all UniRef90 lengths in nucleotides.
    with gzip.open('/data1/gdouglas/db/UniRef/uniref90.fasta.nucleotide_length.tsv.gz', mode='rt') as uniref_lengths_in:
        uniref_lengths_in.readline()
        for uniref_length_line in uniref_lengths_in:
            uniref_length_line = uniref_length_line.rstrip()
            uniref_length_split = uniref_length_line.split('\t')
            uniref90_lengths[uniref_length_split[0]] = int(uniref_length_split[1])

    species = read_ids('/data1/gdouglas/projects/hgt_fragments/mapfiles/species.txt')

    accessions = set()
    for sp in species:
        sp_acc = read_ids('/data1/gdouglas/projects/hgt_fragments/genome_accessions/' + sp + '_accessions.txt')
        for acc in sp_acc:
            accessions.add(acc)

    intergenic_pseudo = set()
    pseudo_sizes = defaultdict(int)
    old2new_ids = dict()
    new2old_ids = dict()
    pseudo_blast_hits = defaultdict(list)
    pseudo_aligned_lengths = defaultdict(list)

    for acc in accessions:

        pseudo_fasta = '/data1/gdouglas/projects/hgt_fragments/pseudogenes/pseudofinder_output/pseudos.fasta/' + acc + '_pseudos.fasta.gz'
        pseudo_gff = '/data1/gdouglas/projects/hgt_fragments/pseudogenes/pseudofinder_output/pseudos.gff/' + acc + '_pseudos.gff.gz'
        intergenic_pseudo_blastx = '/data1/gdouglas/projects/hgt_fragments/pseudogenes/pseudofinder_output/intergenic.fasta.blastX_output/' + acc + '_intergenic.fasta.blastX_output.tsv.gz'

        with gzip.open(pseudo_fasta, mode="rt") as pseudo_fasta_in:
            for pseudo_fasta_line in pseudo_fasta_in:
                if pseudo_fasta_line[0] == '>':
                    last_id = pseudo_fasta_line.split()[0][1:]
                else:
                    pseudo_sizes[last_id] = pseudo_sizes[last_id] + len(pseudo_fasta_line.rstrip())

        with gzip.open(pseudo_gff, mode="rt") as pseudo_gff_in:

            for pseudo_gff_line in pseudo_gff_in:

                if pseudo_gff_line[0] == "#":
                    continue

                pseudo_gff_line = pseudo_gff_line.rstrip()

                pseudo_gff_linesplit = pseudo_gff_line.split('\t')

                if 'Reason(s):Intergenic region with' not in pseudo_gff_linesplit[8]:
                    continue

                if ';old_locus_tag=' not in pseudo_gff_linesplit[8]:
                    sys.exit('Stopping - no old locus tag!')

                pseudo_id_split = pseudo_gff_linesplit[8].split(';')

                old_locus_tag = pseudo_id_split[-1].replace('old_locus_tag=', '')
                pseudo_id = pseudo_id_split[-2].replace('locus_tag=', '')

                old2new_ids[old_locus_tag] = pseudo_id
                new2old_ids[pseudo_id] = old_locus_tag
                intergenic_pseudo.add(pseudo_id)

        # Read through BLASTX-formatted output.
        # These are the columns:
        # qseqid sseqid pident slen mismatch gapopen qstart qend sstart send evalue bitscore stitle
        with gzip.open(intergenic_pseudo_blastx, mode="rt") as intergenic_pseudo_blastx_in:

            for intergenic_blastx_line in intergenic_pseudo_blastx_in:

                line_split = intergenic_blastx_line.split('\t')

                if line_split[0] not in old2new_ids.keys():
                    continue

                qend = int(line_split[7])
                qstart = int(line_split[6])
                if qend >= qstart:
                    pseudo_aligned_lengths[line_split[0]] = pseudo_aligned_lengths[line_split[0]] + [qend - qstart + 1]
                else:
                    pseudo_aligned_lengths[line_split[0]] = pseudo_aligned_lengths[line_split[0]] + [qstart - qend + 1]

                pseudo_blast_hits[line_split[0]] = pseudo_blast_hits[line_split[0]] + [uniref90_lengths[line_split[1]]]

    print('\t'.join(['pseudo_id', 'intergenic_id', 'pseudo_size', 'num_matches', 'mean_aligned_length',
                     'sd_aligned_length', 'mean_matched.ref_length', 'sd_matched.ref_length']))

    for pseudo_id in intergenic_pseudo:

        pseudo_size = pseudo_sizes[pseudo_id]

        old_id = new2old_ids[pseudo_id]

        num_matches = len(pseudo_aligned_lengths[old_id])

        if num_matches == 0:
            sys.exit('For some reason this pseudogene has no matches, but it should! ' + pseudo_id + ' ' + old_id)

        mean_aligned_length = numpy.mean(pseudo_aligned_lengths[old_id])
        sd_aligned_length = numpy.std(pseudo_aligned_lengths[old_id])

        mean_matched_ref_length = numpy.mean(pseudo_blast_hits[old_id])
        sd_matched_ref_length = numpy.std(pseudo_blast_hits[old_id])

        print('\t'.join([pseudo_id,
                         old_id,
                         str(pseudo_size),
                         str(num_matches),
                         str(float("{:.2f}".format(mean_aligned_length))),
                         str(float("{:.2f}".format(sd_aligned_length))),
                         str(float("{:.2f}".format(mean_matched_ref_length))),
                         str(float("{:.2f}".format(sd_matched_ref_length)))]))


if __name__ == '__main__':
    main()
