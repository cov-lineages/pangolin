#!/usr/bin/env python3

import argparse
import collections
from Bio import SeqIO
import os
import sys 
import csv
import lineage_classes

def parse_args():
    parser = argparse.ArgumentParser(description='SNP based classification.')

    parser.add_argument("-snps", action="store", type=str, dest="snps")
    parser.add_argument("-o","--outfile",action="store", type=str, dest="outfile")
    return parser.parse_args()

def traverse(lineages,query,result=("",0,0),lineage="root"):
    query_snps = query.snps
    if query_snps == []:
        query.update_lineage(("A",100,0))
    else:
        for child in lineages[lineage].children:

            lineage_snps = child.snps
            overlap = set(lineage_snps)&set(query_snps)

            if overlap:
                identity = len(overlap)/len(lineage_snps)
                if identity >= result[1]:
                    if len(overlap) > result[2]:
                        prop_included = round(len(overlap)/len(query_snps), 2)
                        result = (child.name, identity, len(overlap))
                        query.update_lineage(result)

            traverse(lineages,query,result,child.name)

def get_ancestors(lineages):
    for l in lineages:
        lineage = lineages[lineage]
        lineage_list = lineage.name.split(".")
        for i in range(len(lineage_list)):
            
            ancestor_lineage = ".".join(lineage_list[:i+1])
            if ancestor_lineage not in lineage.ancestors:
                lineage.ancestors.append(ancestor_lineage)

def add_snps(lineages,snp_dict,lineage="root"):
    for child in lineages[lineage].children:
        for ancestor in child.ancestors:
            snps = snp_dict[ancestor.name]
            for snp in snps:
                child.add_snp(snp)

        recursively_add_snps(lineages,snp_dict,child.name)

def create_lineage_tree(defining_snps):
    lineages = {}
    snp_dict = {}
    with open(defining_snps,"r") as f:
        for l in f:
            l = l.rstrip("\n")
            if not l.startswith("lineage"):
                lineage,snps=l.split(",")
                if lineage not in lineages:
                    this_lineage= Lineage(lineage, lineages)
                    snp_dict[lineage] = snps.split(";")
                    lineages[lineage]= this_lineage
    
    recursively_add_snps(lineages,snp_dict)

    return lineages

def try_to_classify_by_snps():

    args = parse_args()
    
    lineages= create_lineage_tree(args.defining_snps)
    
    queries = []

    with open(args.query_snps,"r") as f:
        for l in f:
            l = l.rstrip("\n")
            name,snp_string= l.split(',')
            snps = snp_string.split(";")

            new_query = Query(name,snps)
            traverse(lineages,new_query)
            queries.append(new_query)

    with open(args.outfile, "w") as fw:
        fw.write("name,lineage,identity,query_snps,num_query_snps_included\n")
        for query in queries:
            query_snp_string = ';'.join(query.snps)
            fw.write(f"{query.name},{query.lineage},{query.identity},{query_snp_string},{query.prop_snps_included}\n")

if __name__ == '__main__':

    try_to_classify_by_snps()