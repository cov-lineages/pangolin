#!/usr/bin/env python3

import argparse
import collections
from Bio import SeqIO
import os
import sys 
import csv

class Query:
    def __init__(self, name, snps):
        self.name = name
        self.snps = snps
        self.lineage = ""
        self.identity = 0
        self.prop_snps_included = 0
        
    def update_lineage(self, lineage):
        self.lineage = lineage[0]
        self.identity = lineage[1]
        self.prop_snps_included = lineage[2]

class Lineage:
    def __init__(self, name, lineages):
        self.name = name
        
        self.snps = []
        self.ancestors = []
        if self.name != "root":
            self.parent = self.assign_parent(lineages)
            self.get_ancestors(lineages)
        self.children = []

    def get_children_names(self):
        print([i.name for i in self.children])
    
    def assign_parent(self, lineages):
        
        name_list = self.name.split(".")
        
        parent_lineage = ""
        if len(name_list) == 1:
            if self.name == "A":
                parent_lineage = "root"
            elif self.name == "B":
                parent_lineage = "A"
        else:
            parent_lineage = '.'.join(self.name.split(".")[:-1])
            
        if parent_lineage in lineages:
            siblings = lineages[parent_lineage].children
            siblings.append(self)
            return lineages[parent_lineage]
            
        else:
            new_parent = Lineage(parent_lineage,lineages)
            siblings = new_parent.children
            siblings.append(self)
            lineages[new_parent.name]=new_parent
            return new_parent
        
    def get_parent_snps(self,snp):
        snps = []
        parent= self.parent
        for snp in parent.snps:
            snps.append(snp)

        if parent.name != "root":
            parent.get_parent_snps(snp)
        return list(set(snps))
    
    def add_snp(self, snp):
        if snp not in self.snps:
            parents_snps = self.get_parent_snps(snp)
            if snp not in parents_snps:
                self.snps.append(snp)

def parse_args():
    parser = argparse.ArgumentParser(description='SNP based classification.')

    parser.add_argument("-snps","--defining-snps", action="store", type=str, dest="defining_snps")
    parser.add_argument("-q","--query-snps", action="store", type=str, dest="query_snps", help="A fasta file containing the query sequences")
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

def recursively_add_snps(lineages,snp_dict,lineage="root"):
    for child in lineages[lineage].children:
        
        snps = snp_dict[child.name]
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