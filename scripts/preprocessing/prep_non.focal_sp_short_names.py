#!/usr/bin/python3

import argparse
import sys
import os

sp_map = dict()
short = set()

ids_file = '/home/gdouglas/projects/def-shapiro-ab/gdouglas/projects/pseudogenes/mapfiles/focal_and_non.focal_species.txt'

with open(ids_file, 'r') as ids_fh:
    for line in ids_fh:
        full = line.rstrip()
        name_split = full.split('_')
        potential_short_basis = []
        for i in range(len(name_split)):
            potential_short_basis.append(name_split[i][0])

        potential_short = '.'.join(potential_short_basis)
        element_to_expand = 1
        str_index_to_add = 1

        while potential_short in short:
            while len(name_split[element_to_expand]) <= str_index_to_add:
                element_to_expand += 1
                str_index_to_add = 1

            potential_short_basis[element_to_expand] = potential_short_basis[element_to_expand] + name_split[element_to_expand][str_index_to_add]
            str_index_to_add += 1
            potential_short = '.'.join(potential_short_basis)

        short.add(potential_short)
        sp_map[potential_short] = full

for k, v in sp_map.items():
    print(v + '\t' + k)