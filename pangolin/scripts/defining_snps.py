import argparse
import collections
from Bio import AlignIO
import os
cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Find lineage defining snps.')

    parser.add_argument("-a", action="store", type=str, dest="a")

    parser.add_argument("-o", action="store", type=str, dest="o")

    return parser.parse_args()

def lineage_dict(alignment_file):
    aln = AlignIO.read(alignment_file, "fasta")

    lineages_dict = collections.defaultdict(list)
    for record in aln:
        if record.id == "outgroup_A":
            lineages_dict["reference"].append(record)
        else:
            lineage = record.id.split("_")[1]
            lineages_dict[lineage].append(record)
    return lineages_dict

def find_snps(ref,member):
    snps = []
    index = 0 
    for i in range(len(ref)):
        if ref[i]!= '-':
            index +=1
            
        col = [ref[i],member[i]]
        if len(set(col))>1:
            if col[1].lower() == 'n':
                pass
            else:
                snp = f"{index}{col[0].upper()}{col[1].upper()}"
                snps.append(snp)
    return snps
                
def get_all_snps_in_lineages(lineage_dict):
    reference = lineage_dict["reference"][0]
    lineage_defining_snps = {}
    for lineage in sorted(lineage_dict):
        if lineage != "reference":
            lineage_snps = []
            for member in lineage_dict[lineage]:
                snps = find_snps(reference.seq, member.seq)
                lineage_snps.append(snps)
            lineage_set = list(set.intersection(*[set(x) for x in lineage_snps]))
            lineage_list = sorted(lineage_set, key = lambda x : int(x[:-2]))
            snp_str = ";".join(lineage_list)

            lineage_defining_snps[lineage] = snp_str
            
            if lineage_list != []:
                print(f"Lineage defining snps for {lineage}")
                for i in lineage_list:
                    print("\t",i)
            else:
                print(f"No lineage defining snps found for {lineage}")

    return lineage_defining_snps

def read_alignment_and_get_snps():

    args = parse_args()

    alignment_file = os.path.join(cwd, args.a)
    if not os.path.exists(alignment_file):
        sys.stderr.write('Error: cannot find alignment file at {}\n'.format(alignment_file))
        sys.exit(-1)
    else:
        print(f"Reading in alignment file {alignment_file}.")

    ldict = lineage_dict(alignment_file)
    sdict = get_all_snps_in_lineages(ldict)

    with open(args.o, "w") as fw:
        fw.write("lineage,defining_snps\n")
        for lineage in sdict:
            snps = sdict[lineage]
            fw.write(f"{lineage},{snps}\n")

if __name__ == '__main__':

    read_alignment_and_get_snps()