#!/usr/bin/env python3

import argparse
import collections
from Bio import AlignIO
import os
cwd = os.getcwd()

"""
Current output:
--snps
taxon,snps
taxon1,2897GT;30000TA

"""


def parse_args():
    parser = argparse.ArgumentParser(description='Find all snps.')

    parser.add_argument("-a", action="store", type=str, dest="a")

    parser.add_argument("-o", action="store", type=str, dest="snps")
    return parser.parse_args()

def find_snps(ref,member):
    """Identifies unambiguous snps between two sequences 
    and returns them as a list, using position in the ref seq (i.e. no gaps in ref)"""
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

def get_reference(alignment,reference_id="Wuhan/WH04/2020"):
    """get reference seq record"""
    for record in alignment:
        if record.id == reference_id:
            return record

def pcent_done(c, total):
    return round((c*100)/total, 2)

def write_snps(alignment, reference,outfile):
    """add list of snps relative to ref as an annotation to the seq record"""
    c = 0
    total = len(alignment)
    for record in alignment:
        c +=1 
        if c%500==0:
            print(pcent_done(c, total), '%')

        snps = find_snps(reference.seq, record.seq)
        snp_string =snp_list_to_snp_string(snps)
        outfile.write(f"{record.id},{snp_string}\n")

    print(total, "records annotated")

def get_ids_in_list_of_records(records):
    """return ids in a set of seq records"""
    ids = []
    for record in records:
        ids.append(record.id)
    return ids

def snp_list_to_snp_string(snp_list):
    """turn a snp list into a `;`-separated string of snps that are sorted by 
    position in the genome"""
    snp_string = ";".join(sorted(snp_list, key = lambda x : int(x[:-2])))
    return snp_string

def get_all_snps(alignment_file,outfile):
    """ this is the main worker function of this script. 
    ultimately it returns a list of singleton snps to_mask
    and a list of lineage_defining_snps per lineage to write to a file

    1. reads in the alignment
    2. identifies the reference sequence
    3. adds in some useful things to the seq record (n, snps, snp_string)

    4. structures alignment records by lineage
    5. for each lineage
        - count occurences of each snp and id singletons to mask
        - make dict keyed by unique snp combinations with all the associated records as value list
        - get snps to be represented in the tree and potential defining snps (by % cut off args)
        - get best taxa to represent the snps to be represented based on N content
        - if theres too few taxa for each lineage add some more (total of 5)
        - write representatives 
        - get a set of snps that 100 of taxa are in
        - if snp set is empty, pad with snps flagged as potential defining snps by the less strict cut off
    6. return to_mask and lineage_defining_snps
    """
    print("1. Reading in the alignment")
    aln = AlignIO.read(alignment_file, "fasta")
    print("2. Getting the reference:")
    reference = get_reference(aln)
    print(reference.id)

    print("3. Find and write all snps")
    write_snps(aln, reference,outfile)


def read_alignment_and_write_files():

    args = parse_args()

    alignment_file = os.path.join(cwd, args.a)
    if not os.path.exists(alignment_file):
        sys.stderr.write('Error: cannot find alignment file at {}\n'.format(alignment_file))
        sys.exit(-1)
    else:
        print(f"Reading in alignment file {alignment_file}.")

    with open(args.snps,"w") as fw:
        fw.write("taxon,snps\n")
        get_all_snps(alignment_file,fw)


if __name__ == '__main__':

    read_alignment_and_write_files()