#!/usr/bin/env python3

import argparse
import collections
from Bio import AlignIO
import os
import csv
import dendropy

cwd = os.getcwd()


def parse_args():
    parser = argparse.ArgumentParser(description='Find all snps.')

    parser.add_argument("--lineages", action="store", type=str, dest="l")
    parser.add_argument("--metadata",action="store", type=str, dest="metadata")
    parser.add_argument('--polytomies',action="store", type=str, dest="polytomies")


    parser.add_argument("--outfile", action="store", type=str, dest="outfile")

    return parser.parse_args()


def add_phylotype_annotation(alignment,metadata):
    phylotype = {}
    with open(metadata,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            phylotype[row["sequence_name"]]=row["phylotype"]
    for record in alignment:
        record.annotations["phylotype"] = phylotype[record.id]

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


def get_phylotypes(metadata):
    phylotype = {}
    with open(metadata,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            phylotype[row["sequence_name"]]=row["phylotype"]

    return phylotype

def get_polytomy(polytomies):

    polytomy = {}
    with open(polytomies,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            polytomy[row["name"]]=row["polytomy"]
    return polytomy


def get_lineage_dict(lineage_file):

    lineages = collections.defaultdict(list)

    with open(lineage_file,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            tax_dict = {"name":row["taxon"]}
            lineages[row["lineage"]].append(tax_dict)

    return lineages

def find_basal_polytomy(polytomies, metadata, lineage_file, outfile):
    print("Getting phylotypes")
    phylotypes = get_phylotypes(metadata)
    print(f"{len(phylotypes)} read in")
    print("Getting polytomies")
    polytomy_dict = get_polytomy(polytomies)
    print(f"{len(polytomy_dict)} members assigned")
    print("Getting lineage dict")
    lineages = get_lineage_dict(lineage_file)
    print(f"{len(lineages)} lineages loaded")


    for lineage in lineages:
        print("\nLineage", lineage)

        lineage_taxa = lineages[lineage]
        print(len(lineage_taxa), "members")

        within_lineage_phylotypes = []
        for record in lineage_taxa:

            if record["name"] in phylotypes:
                phylotype =phylotypes[record["name"]]
                record["phylotype"] = phylotype
                if phylotype == "":
                    length = 0
                else:
                    length = len(phylotype.split("."))
                within_lineage_phylotypes.append((phylotype, length))
        print("Within lineage phylotypes:")
        shortest_phylotype= sorted(within_lineage_phylotypes, key = lambda x : x[1])[0]
        print(shortest_phylotype)
        basal_phylotype = []
        for record in lineage_taxa:
            if record["name"] in phylotypes:
                if shortest_phylotype[0] == record["phylotype"]:
                    
                    basal_phylotype.append(record["name"])

        basal_polytomy = []
        for record in basal_phylotype:
            if record in polytomy_dict:
                basal_polytomy.append(record)
        print("Basal polytomy", len(basal_polytomy))

        if basal_polytomy == []:
            print("Switching to basal non-polytomy")
            basal_polytomy = basal_phylotype

        for record in basal_polytomy:
            print("Taxon", record)
            outfile.write(f"{lineage},{record}\n")

def read_in_data_get_basal_polytomy():

    args = parse_args()
    polytomies = os.path.join(cwd, args.polytomies)
    if not os.path.exists(polytomies):
        sys.stderr.write('Error: cannot find polytomy file at {}\n'.format(polytomies))
        sys.exit(-1)
    else:
        print("Tree file found at", polytomies)

    metadata_file = os.path.join(cwd, args.metadata)
    if not os.path.exists(metadata_file):
        sys.stderr.write('Error: cannot find metadata file at {}\n'.format(metadata_file))
        sys.exit(-1)
    else:
        print(f"Reading in metadata file {metadata_file}.")

    
    lineage_file = os.path.join(cwd, args.l)
    if not os.path.exists(lineage_file):
        sys.stderr.write('Error: cannot find lineage file at {}\n'.format(lineage_file))
        sys.exit(-1)
    else:
        print(f"Reading in lineage annotations file {lineage_file}.")

    with open(args.outfile, "w") as fw:
        fw.write("lineage,taxon\n")

        find_basal_polytomy(polytomies,metadata_file,lineage_file, fw)


if __name__ == '__main__':

    read_in_data_get_basal_polytomy()