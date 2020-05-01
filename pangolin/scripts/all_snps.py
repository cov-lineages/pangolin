#!/usr/bin/env python3

import argparse
import collections
from Bio import AlignIO
import os
cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Find all snps.')

    parser.add_argument("-a", action="store", type=str, dest="a")
    parser.add_argument("-l", action="store", type=str, dest="l")

    parser.add_argument("--all_snps", action="store", type=str, dest="o")
    parser.add_argument("--defining_snps", action="store", type=str, dest="d")
    parser.add_argument("--to_mask", action="store", type=str, dest="m")
    return parser.parse_args()

def get_lineage_dict(alignment_file, lineage_file):
    """Takes in lineage annotations and an alignment file. 
    Outputs structured dict of seq records per lineage."""
    aln = AlignIO.read(alignment_file, "fasta")
    lineages_dict = {}
    lineages_records = collections.defaultdict(list)

    with open(lineage_file, "r") as f:
        for l in f:
            l = l.rstrip("\n")
            tokens = l.split(",")
            lineages_dict[tokens[0]]=tokens[1]

    for record in aln:
        if record.id == "Wuhan/WH04/2020":
            lineages_records["reference"].append(record)
        else:
            try:
                lineage = lineages_dict[record.id]
                lineages_records[lineage].append(record)
            except:
                print(record.id, "found in alignment but not lineages csv file")

    print("Lineage\t\tNum sequences")
    for lineage in sorted(lineages_records):
        if lineage != "reference":
            print(f"{lineage}\t\t{len(lineages_records[lineage])}")
    
    return lineages_records

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

def get_N_content(seq):
    num_N = str(seq).upper().count("N")
    pcent_N = (num_N*100)/len(seq)
    return pcent_N

def get_all_snps(alignment_file,lineage_file,outfile,defining_file,fmask):
    fm = open(fmask,"w")
    fd = open(defining_file,"w")
    lineages_dict = get_lineage_dict(alignment_file, lineage_file)
    reference = lineages_dict["reference"][0]
    
    for lineage in sorted(lineages_dict):
        if lineage != "reference":
            seq_dict = {}
            n_dict = {}
            lineage_snps = collections.defaultdict(list)
            member_snps = {}
            snp_counter = collections.defaultdict(list)

            for member in lineages_dict[lineage]:
                seq_dict[member.id] = member.seq
                snps = find_snps(reference.seq, member.seq)
                for i in snps:
                    snp_counter[i].append(member.id)

                sorted_snps = ";".join(sorted(snps, key = lambda x : int(x[:-2])))
                member_snps[member.id]= snps
                
                pcent_N = get_N_content(member.seq)
                n_dict[member.id]=pcent_N
                lineage_snps[sorted_snps].append((member.id,pcent_N))

            potential_lineage_defining = []
            flagged = []

            for snp in snp_counter:
                if len(snp_counter[snp]) ==1:
                    print("singleton",snp)
                    fm.write(f"{lineage},{snp}\n")
                inclusion_pcent = (100*len(snp_counter[snp])/len(lineages_dict[lineage]))
                if inclusion_pcent > 90:
                    potential_lineage_defining.append(snp)
                if inclusion_pcent > 10:
                    flagged.append(snp)

            represented= {}
            for snp_set in lineage_snps:
                snps = snp_set.split(";")
                
                for snp in snps:
                    if snp in flagged and not snp in represented:
                        member = sorted(lineage_snps[snp_set], key = lambda x : int(x[1]))[0]
                        represented[snp] = member[0]
                        for i in member_snps[member[0]]:
                            represented[i] = member[0]

            representative_taxa = []
            for snp in represented:
                representative_taxa.append(represented[snp])
            representative_taxa=list(set(representative_taxa))

            if len(representative_taxa)<5:
                take_best_of_rest = []
                members = lineages_dict[lineage]
                for member in members:
                    if member.id not in representative_taxa:
                        take_best_of_rest.append((member.id, n_dict[member.id]))
                take_best_of_rest = sorted(take_best_of_rest, key = lambda x : x[1])
                
                for item in take_best_of_rest:
                    if len(representative_taxa) < 5:
                        if item[0] not in representative_taxa:
                            representative_taxa.append(item[0])
            
            print(f"{lineage}: {len(representative_taxa)} representative seqs")

            for i in representative_taxa:
                snps = member_snps[i]
                snp_string = ";".join(sorted(snps, key = lambda x : int(x[:-2])))
                outfile.write(f"{lineage},{i},{snp_string}\n")
            

            lineage_set = list(set.intersection(*[set(x.split(";")) for x in lineage_snps]))
            print(f"{lineage} Lineage set:",lineage_set)
            if lineage_set == []:
                for i in potential_lineage_defining:
                    lineage_set.append(i)
            lineage_str = ";".join(sorted(lineage_set, key = lambda x : int(x[:-2])))
            fd.write(f"{lineage},{lineage_str}\n")
    fm.close()
    fd.close()

def read_alignment_and_get_snps():

    args = parse_args()

    alignment_file = os.path.join(cwd, args.a)
    if not os.path.exists(alignment_file):
        sys.stderr.write('Error: cannot find alignment file at {}\n'.format(alignment_file))
        sys.exit(-1)
    else:
        print(f"Reading in alignment file {alignment_file}.")

    lineage_file = os.path.join(cwd, args.l)
    if not os.path.exists(lineage_file):
        sys.stderr.write('Error: cannot find lineage file at {}\n'.format(lineage_file))
        sys.exit(-1)
    else:
        print(f"Reading in lineage annotations file {lineage_file}.")

    with open(args.o, "w") as fw:
        fw.write("lineage,name,snps\n")
        get_all_snps(alignment_file,lineage_file,fw,args.d,args.m)

if __name__ == '__main__':

    read_alignment_and_get_snps()