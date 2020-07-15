#!/usr/bin/env python3

import argparse
import collections
from Bio import AlignIO
import os
import csv
import dendropy

cwd = os.getcwd()


def parse_args():
    parser = argparse.ArgumentParser(description='Find all polytomies.')

    parser.add_argument('--global-tree',action="store", type=str, dest="global_tree")

    parser.add_argument("--outfile", action="store", type=str, dest="outfile")

    return parser.parse_args()

def process_node(node, polytomies):
    if len(node.child_nodes()) > 0:
        zbl = []
        for child in node.child_nodes():
            length = int(child.edge.length * 29903)
            if (length == 0 and len(child.child_nodes()) == 0):
                zbl.append(child.taxon.label)
        if (len(zbl) > 1):
            for taxon in zbl:
                polytomies[taxon] = zbl
        for child in node.child_nodes():
            process_node(child, polytomies)


def get_polytomy(global_tree):

    tree = dendropy.Tree.get(file=open(global_tree, "r"), schema="nexus")

    polytomies = {}
    process_node(tree.seed_node, polytomies)
    
    return polytomies


def read_in_data_get_polytomy():

    args = parse_args()
    global_tree_file = os.path.join(cwd, args.global_tree)
    if not os.path.exists(global_tree_file):
        sys.stderr.write('Error: cannot find global tree file at {}\n'.format(global_tree_file))
        sys.exit(-1)
    else:
        print("Tree file found at", global_tree_file)

    fb = open(args.outfile, "w")
    fb.write("name,polytomy\n")

    polytomies = get_polytomy(global_tree_file)
    for tax in polytomies:

        polytomies_str = ";".join(polytomies[tax])
        fb.write(f"{tax},{polytomies_str}\n")
    fb.close()


if __name__ == '__main__':

    read_in_data_get_polytomy()