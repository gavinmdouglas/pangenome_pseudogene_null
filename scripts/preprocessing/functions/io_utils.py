#!/usr/bin/python3

import gzip
import sys
import textwrap
from collections import OrderedDict

def read_ids(filename, split_line = False):
    '''Read ids from a file into a list (one id per line)'''
    ids = list()
    with open(filename, 'r') as id_file:
        for line in id_file:
            if split_line:
                ids.append(line.split()[0])
            else:
                ids.append(line.rstrip())

    return(ids)


def read_fasta(filename, cut_header=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value.'''

    # Intitialize empty dict.
    seq = {}

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:

        line = line.rstrip()

        if len(line) == 0:
            continue

        # If header-line then split by whitespace, take the first element,
        # and define the sequence name as everything after the ">".
        if line[0] == ">":

            if cut_header:
                name = line.split()[0][1:]
            else:
                name = line[1:]

            name = name.rstrip("\r\n")

            # Make sure that sequence id is not already in dictionary.
            if name in seq:
                sys.stderr("Stopping due to duplicated id in file: " + name)

            # Intitialize empty sequence with this id.
            seq[name] = ""

        else:
            # Remove line terminator/newline characters.
            line = line.rstrip("\r\n")

            # Add sequence to dictionary.
            seq[name] += line

    fasta_in.close()

    return seq


def read_ordered_fasta(filename, cut_header=False):
    '''Read in FASTA file (gzipped or not) and return dictionary with each
    independent sequence id as a key and the corresponding sequence string as
    the value. Note that the type will be OrderedDict to maintain the order
    of sequences in the file.'''

    # Intitialize empty dict.
    seq = OrderedDict()

    # Intitialize undefined str variable to contain the most recently parsed
    # header name.
    name = None

    # Read in FASTA line-by-line.
    if filename[-3:] == ".gz":
        fasta_in = gzip.open(filename, "rt")
    else:
        fasta_in = open(filename, "r")

    for line in fasta_in:

        line = line.rstrip()

        if len(line) == 0:
            continue

        # If header-line then split by whitespace, take the first element,
        # and define the sequence name as everything after the ">".
        if line[0] == ">":

            if cut_header:
                name = line.split()[0][1:]
            else:
                name = line[1:]

            name = name.rstrip("\r\n")

            # Make sure that sequence id is not already in dictionary.
            if name in seq:
                sys.stderr("Stopping due to duplicated id in file: " + name)

            # Intitialize empty sequence with this id.
            seq[name] = ""

        else:
            # Remove line terminator/newline characters.
            line = line.rstrip("\r\n")

            # Add sequence to dictionary.
            seq[name] += line

    fasta_in.close()

    return seq


def write_fasta(seq, outfile):
    out_fasta = open(outfile, "w")

    # Look through sequence ids (sorted alphabetically so output file is
    # reproducible).
    for s in sorted(seq.keys()):
        out_fasta.write(">" + s + "\n")
        out_fasta.write(textwrap.fill(seq[s], width=70) + "\n")

    out_fasta.close()


def read_fastq_headers(filepath):
        
        # Capture every 4th line of FASTQ (after first line)
        lineno = 4
        header_lines = []

        with open(filepath, "r") as fastq_in:
            
            for line in fastq_in:
                
                if lineno == 4:
                    header_lines.append(line.rstrip())
                    lineno = 1
                
                else:
                    lineno += 1

        return(header_lines)

