#!/usr/bin/python3

from collections import defaultdict

### Get map of cd-hit-est clusters to representative gene ids.

def main():

    with open('/data1/gdouglas/projects/accessory_vs_pseudogene/cdhit_workflow_intact_and_intergenic/cdhit_filt_intergenic_pseudogenes_and_intact_c0.95_aL0.90_s0.90.clstr', 'r') as cdhitest_out:
        for line in cdhitest_out:
            line_split = line.split()

            if len(line_split) == 0:
                continue

            if line_split[0][0] == ">":
                current_cluster = 'c' + line_split[1]

            if line_split[-1] == "*":
                gene = line_split[2][1:]
                gene = gene.replace(".", "")
                print(current_cluster + "\t" + gene)


if __name__ == '__main__':
    main()
