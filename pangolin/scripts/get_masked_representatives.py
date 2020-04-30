#!/usr/bin/env python3

import argparse
import collections
from Bio import AlignIO
from Bio import SeqIO
import os
import sys 
import csv

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Extract representative seqs with masked singleton snps.')

    parser.add_argument("-r", action="store", type=str, dest="r")
    parser.add_argument("-m", action="store", type=str, dest="m")
    parser.add_argument("-l", action="store", type=str, dest="l")
    parser.add_argument("-a", action="store", type=str, dest="a")
    parser.add_argument("-o", action="store", type=str, dest="o")
    parser.add_argument("--metadata", action="store", type=str, dest="mi")
    parser.add_argument("--metadata-out", action="store", type=str, dest="mo")
    return parser.parse_args()

def find_snps(ref,member):
    """Identifies unambiguous snps between two sequences 
    and returns them as a list"""
    snps = []
    index = 0 
    for i in range(len(ref)):
        if ref[i]!= '-':
            index +=1
            
        col = [ref[i],member[i]]
        if len(set(col))>1:
            if not col[1].lower() in ["a","g","t","c","-"]:
                pass
            else:
                snp = f"{index}{col[0].upper()}{col[1].upper()}"
                snps.append(snp)
    return snps

def mask_snp(ref,member,to_mask_snp):
    """Identifies unambiguous snps between two sequences 
    and returns them as a list"""
    new_str = ""
    index = 0 
    for i in range(len(ref)):
        if ref[i]!= '-':
            index +=1
            
        col = [ref[i],member[i]]
        if len(set(col))>1:
            if not col[1].lower() in ["a","g","t","c","-"]:
                pass
                new_str +=col[1]
            else:
                snp = f"{index}{col[0].upper()}{col[1].upper()}"
                list_seq = list(member)
                if snp == to_mask_snp:
                    list_seq[i] = "N"
                    new_str += "N"
                    member = "".join(list_seq)
                else:
                    new_str += col[1]
        else:
            new_str += col[1]
    return new_str

def make_lineage_dict(lineage):
    lineages = {}
    with open(lineage,"r") as f:
        for l in f:
            tokens = l.rstrip().split(",")
            lineages[tokens[0]]=tokens[1]
    return lineages
        
def make_rep_dict(r):
    reps = {}
    with open(r,"r") as f:
        for l in f:
            tokens = l.rstrip().split(",")
            reps[tokens[1]]=tokens[0]
    return reps
        
def make_mask_dict(m):
    to_mask = collections.defaultdict(list)
    with open(m,"r") as f:
        for l in f:
            tokens = l.rstrip().split(",")
            to_mask[tokens[0]].append(tokens[1])
    return to_mask

def get_reference(fasta):
    reference = ""
    for record in SeqIO.parse(fasta,"fasta"):
        if "WH04" in record.id:
            reference = record
    return reference

def extract_representatives_and_do_the_masking_thing():

    args = parse_args()

    rep_file = os.path.join(cwd, args.r)
    if not os.path.exists(rep_file):
        sys.stderr.write('Error: cannot find rep file at {}\n'.format(rep_file))
        sys.exit(-1)
    else:
        print(f"Reading in rep file {rep_file}.")
        reps = make_rep_dict(rep_file)
        print("Number of representative seqs:",len(reps))

    mask_file = os.path.join(cwd, args.m)
    if not os.path.exists(mask_file):
        sys.stderr.write('Error: cannot find mask file at {}\n'.format(mask_file))
        sys.exit(-1)
    else:
        print(f"Reading in mask file {mask_file}.")
        to_mask = make_mask_dict(mask_file)

    lineage_file = os.path.join(cwd, args.l)
    if not os.path.exists(lineage_file):
        sys.stderr.write('Error: cannot find lineage file at {}\n'.format(lineage_file))
        sys.exit(-1)
    else:
        print(f"Reading in lineage file {lineage_file}.")
        lineage_dict = make_lineage_dict(lineage_file)

    metadata = os.path.join(cwd, args.mi)
    if not os.path.exists(metadata):
        sys.stderr.write('Error: cannot find metadata file at {}\n'.format(metadata))
        sys.exit(-1)
    else:
        print(f"Reading in metadata file {metadata}.")
        try:
            with open(metadata,newline="") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    sequence_name = row["sequence_name"]
                    covv_accession_id = row["covv_accession_id"]
                    edin_admin_0 = row["edin_admin_0"]
                    edin_travel = row["edin_travel"]
                    covv_collection_date = row["covv_collection_date"]
                    edin_epi_week = row["edin_epi_week"]
        except:
            sys.stderr.write('Error: unexpected headers on {}\n'.format(metadata))
            sys.exit(-1)


    alignment_file = os.path.join(cwd, args.a)
    reference = ''
    if not os.path.exists(alignment_file):
        sys.stderr.write('Error: cannot find alignment file at {}\n'.format(alignment_file))
        sys.exit(-1)
    else:
        print(f"Reading in alignment file {alignment_file}.")
        reference = get_reference(alignment_file)
        print(f"Using {reference.id} as reference")

    fw = open(args.o,"w")
    fw.write(f">{reference.id}|A\n{reference.seq}\n")
    for record in AlignIO.read(alignment_file,"fasta"):
        if record.id != reference.id and record.id in reps:
            snps = find_snps(reference.seq,record.seq)
            lineage = reps[record.id]
            
            if lineage in to_mask:
                new_seq = record.seq
                for snp in snps:
                    if snp in to_mask[lineage]:
                        
                        
                        new_seq = mask_snp(reference.seq, new_seq, snp)
                    else:
                        pass
                        
                new_snps = find_snps(reference.seq,new_seq)

                fw.write(f">{record.id}|{lineage}\n{new_seq}\n")
                
            else:
                fw.write(f">{record.id}|{lineage}\n{record.seq}\n")

    fw.close()

    fm = open(args.mo, "w")
    header = f"name,GISAID ID,country,travel history,sample date,epiweek,lineage,representative\n"
    fm.write(header)
    c,cin,r = 0,0,0
    with open(metadata,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sequence_name = row["sequence_name"]
            print(sequence_name)
            covv_accession_id = row["covv_accession_id"]
            edin_admin_0 = row["edin_admin_0"]
            edin_travel = row["edin_travel"]
            covv_collection_date = row["covv_collection_date"]
            edin_epi_week = row["edin_epi_week"]

            if sequence_name in lineage_dict:
                rep = 0
                lineage = lineage_dict[sequence_name]
                if sequence_name in reps:
                    rep = 1
                    r +=1
                new_l = f"{sequence_name},{covv_accession_id},{edin_admin_0},{edin_travel},{covv_collection_date},{edin_epi_week},{lineage},{rep}\n"
                cin +=1
                fm.write(new_l)
            else:
                c+=1

    print("Filtered out:",c)
    print("Included in the metadata:", cin)
    print("Number of representatives indicated:",r)
    fm.close()

if __name__ == '__main__':

    extract_representatives_and_do_the_masking_thing()