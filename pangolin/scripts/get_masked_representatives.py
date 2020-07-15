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

    parser.add_argument("--representatives", action="store", type=str, dest="representatives")
    parser.add_argument("--to-mask", action="store", type=str, dest="to_mask")
    parser.add_argument("-l", action="store", type=str, dest="l")
    parser.add_argument("-a", action="store", type=str, dest="a")
    parser.add_argument("--representative-seqs-out", action="store", type=str, dest="representatives_out")
    parser.add_argument("--metadata", action="store", type=str, dest="metadata_in")
    parser.add_argument("--metadata-out", action="store", type=str, dest="metadata_out")
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
                new_str +=col[1]
            else:
                snp = f"{index}{col[0].upper()}{col[1].upper()}"
                if snp == to_mask_snp:
                    new_str += "?"
                else:
                    new_str += col[1]
        else:
            new_str += col[1]
    return new_str

def make_lineage_dict(lineage):
    """return lineage dict lineage_dict[name]=lineage"""
    lineages = {}
    with open(lineage,"r") as f:
        for l in f:
            tokens = l.rstrip().split(",")
            lineages[tokens[0]]=tokens[1]
    return lineages
        
def make_rep_dict(r):
    """
    --representative-seqs
    lineage,name
    B,WH0X/Taxon/Name
    """
    reps = {}
    with open(r,"r") as f:
        for l in f:
            tokens = l.rstrip().split(",")
            reps[tokens[1]]=tokens[0]
    return reps
        
def make_mask_dict(m):
    """
    return mask dict mask_dict[lineage]=snp
    input:
    --to_mask
    lineage,snp,taxon
    B,2897GT,WH0X/Taxon/Name
    """
    to_mask = collections.defaultdict(list)
    with open(m,"r") as f:
        for l in f:
            tokens = l.rstrip().split(",")
            to_mask[tokens[0]].append(tokens[1])
    return to_mask

def get_reference(fasta):
    """return reference seq record """
    reference = ""
    for record in SeqIO.parse(fasta,"fasta"):
        if "Wuhan/WH04/2020" == record.id:
            reference = record
    return reference

def make_masked_representative_fasta(alignment_file, reps, reference, fw, to_mask):
    aln = AlignIO.read(alignment_file,"fasta")
    print("SNPs premask\t\tSNPs to mask\t\tSNPs postmask")
    for record in aln:
        if record.id != reference.id and record.id in reps:
            snps = find_snps(reference.seq,record.seq)
            lineage = reps[record.id]
            snp_mask_count = 0
            if lineage in to_mask:
                
                new_seq = record.seq
                for snp in snps:
                    if snp in to_mask[lineage]:
                        snp_mask_count +=1
                        new_seq = mask_snp(reference.seq, new_seq, snp)
                    else:
                        pass
                
                new_snps = find_snps(reference.seq,new_seq)
                print(f"{len(snps)}\t\t{snp_mask_count}\t\t{len(new_snps)}")

                fw.write(f">{record.id}|{lineage}\n{new_seq}\n")
            else:
                fw.write(f">{record.id}|{lineage}\n{record.seq}\n")

def make_metadata_out(metadata,lineage_dict,reps,metadata_out_file):
    c,cin,r = 0,0,0
    with open(metadata,newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            sequence_name = row["sequence_name"]
            gisaid_id = row["covv_accession_id"]
            country = row["country"]
            travel_history = row["travel_history"]
            sample_date = row["sample_date"]
            epi_week = row["epi_week"]

            if sequence_name in lineage_dict:
                rep = 0
                lineage = lineage_dict[sequence_name]
                if sequence_name in reps:
                    rep = 1
                    r +=1
                new_l = f"{gisaid_id},{sequence_name},{country},{travel_history},{sample_date},{epi_week},{lineage},{rep}\n"
                cin +=1
                metadata_out_file.write(new_l)
            else:
                c+=1

    print("Filtered out:",c)
    print("Included in the metadata:", cin)
    print("Number of representatives indicated:",r)

def extract_representatives_and_do_the_masking_thing():
    """
    input:
    --to_mask
    lineage,snp,taxon
    B,2897GT,WH0X/Taxon/Name

    --defining-snps
    lineage,defining_snps
    B,2897GT;30000TA

    --representative-seqs
    lineage,name
    B,WH0X/Taxon/Name
    """
    args = parse_args()

    rep_file = os.path.join(cwd, args.representatives)
    if not os.path.exists(rep_file):
        sys.stderr.write('Error: cannot find rep file at {}\n'.format(rep_file))
        sys.exit(-1)
    else:
        print(f"Reading in rep file {rep_file}.")
        reps = make_rep_dict(rep_file)
        print("Number of representative seqs:",len(reps))

    mask_file = os.path.join(cwd, args.to_mask)
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

    metadata = os.path.join(cwd, args.metadata_in)
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
                    edin_admin_0 = row["country"]
                    travel_history = row["travel_history"]
                    sample_date = row["sample_date"]
                    epi_week = row["epi_week"]
        except:
            sys.stderr.write("Error: unexpected headers in {}\n. \
                            Expected header names:\n\
                            sequence_name,\
                            country,travel_history,sample_date,\
                            epi_week\n".format(metadata))
            sys.exit(-1)


    alignment_file = os.path.join(cwd, args.a)
    if not os.path.exists(alignment_file):
        sys.stderr.write('Error: cannot find alignment file at {}\n'.format(alignment_file))
        sys.exit(-1)
    else:
        print(f"Reading in alignment file {alignment_file}")

    reference = get_reference(alignment_file)
    print(f"Using {reference.id} as reference")

    with open(args.representatives_out,"w") as fw:
        fw.write(f">{reference.id}|A\n{reference.seq}\n")
        make_masked_representative_fasta(alignment_file, reps, reference, fw, to_mask)

    with open(args.metadata_out, "w") as fm:
        header = f"GISAID ID,name,country,travel history,sample date,epiweek,lineage,representative\n"
        fm.write(header)
        make_metadata_out(metadata,lineage_dict,reps,fm)

if __name__ == '__main__':

    extract_representatives_and_do_the_masking_thing()