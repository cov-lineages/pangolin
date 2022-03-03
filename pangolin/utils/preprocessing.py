import os
import sys
from pangolin.utils.log_colours import green,cyan
import hashlib
from Bio import SeqIO
import csv
import gzip


def create_seq_hash(seq_file,hash_map,hashed_seqs):
    """
    iterates through seq file, keeps track of seq hash
    writes record id -> hash map as it goes along
    writes out records as it goes along
    
    returns record count and collapsed record count
    """
    hash_dict = {}
    record_count = 0

    fw2 = open(hashed_seqs, "w")

    with open(hash_map,"w") as fw:
        fw.write("name\thash\n")
        for record in SeqIO.parse(seq_file, "fasta"):
            record_count +=1
            
            seq = str(record.seq).encode()
            hash_object = hashlib.md5(seq)
            hash_str = hash_object.hexdigest()
            fw.write(f"{record.description}\t{hash_str}\n")

            if hash_str not in hash_dict:
                
                fw2.write(f">{hash_str}\n{record.seq}\n")
                hash_dict[hash_str] = 1

    fw2.close()

    return record_count,len(hash_dict)

def designation_assign(designation_cache,hash_map,outfile):
    """
    loads the designation hash file as a dictionary
    runs through query file- checks if in designation dict
    writes to csv the status
    """
    pango_hash = {}
    with open(designation_cache,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            pango_hash[row["seq_hash"]] = row["lineage"]
    
    designated = 0
    with open(outfile,"w") as fw:
        fw.write("hash,designated,lineage\n")
        with open(hash_map, "r") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                seq_hash = row["hash"]
                if seq_hash in pango_hash:
                    designated +=1
                    
                    lineage = pango_hash[seq_hash]
                    fw.write(f"{seq_hash},True,{lineage}\n")
                else:
                    fw.write(f"{seq_hash},False,\n")

    return designated

def seq_qc(in_fasta,pass_qc,qc_status,max_ambiguity):
    """
    loads the aligned, hashed fasta file
    calculates n percentage
    if it's greater than max ambiguity, qc is failed
    writes failed and passed to a file
    """
    fw_pass = open(pass_qc,"w")
    with open(qc_status,"w") as fw:
        fw.write("hash,qc_status,qc_notes\n")

        total_input = 0
        total_pass = 0
    
        for record in SeqIO.parse(in_fasta, "fasta"):
            total_input +=1

            num_N = str(record.seq).upper().count("N")
            prop_N = round((num_N)/len(record.seq), 2)
            if prop_N > max_ambiguity:
                fw.write(f"{record.id},fail,Ambiguous_content:{prop_N}\n")
            else:
                total_pass +=1

                fw.write(f"{record.id},pass,Ambiguous_content:{prop_N}\n")
                fw_pass.write(f">{record.id}\n{record.seq}\n")
    fw_pass.close()

    return total_pass
    

def merge_files(fasta, qc_status, scorpio_report, designated, hash_map, out_merged):
    """
    gets the status of scorpio assignments, designation assignments, seq qc
    as well as the hashing map
    goes through the original fasta file and for every record collates the 
    available info into a csv file with header fields as below:
    """
    header = ["name","hash","lineage","scorpio_constellations","scorpio_mrca_lineage","scorpio_incompatible_lineages","scorpio_support","scorpio_conflict","scorpio_notes","designated","qc_status","qc_notes"]
    with open(out_merged,"w") as fw:
        writer = csv.DictWriter(fw, fieldnames=header,lineterminator="\n")
        writer.writeheader()
        info_dict = {}
        name_dict = {}

        with open(hash_map, "r") as f:
            reader = csv.DictReader(f,delimiter="\t")
            for row in reader:
                info_dict[row["hash"]] = {
                    "hash":row["hash"]
                    }
                name_dict[row["name"]] = row["hash"]
        
        with open(designated, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                info_dict[row["hash"]]["designated"] = row["designated"]
                info_dict[row["hash"]]["lineage"] = row["lineage"]
            
        with open(qc_status, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                info_dict[row["hash"]]["qc_status"] = row["qc_status"]
                info_dict[row["hash"]]["qc_notes"] = row["qc_notes"]
        
        with open(scorpio_report,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                info_dict[row["query"]]["scorpio_constellations"] = row["constellations"]
                info_dict[row["query"]]["scorpio_mrca_lineage"] = row["mrca_lineage"]
                info_dict[row["query"]]["scorpio_incompatible_lineages"] = row["incompatible_lineages"]
                if row["support"]:
                    info_dict[row["query"]]["scorpio_support"] = round(float(row["support"]),2)
                    info_dict[row["query"]]["scorpio_conflict"] = round(float(row["conflict"]),2)
                else:
                    info_dict[row["query"]]["scorpio_support"] = row["support"]
                    info_dict[row["query"]]["scorpio_conflict"] = row["conflict"]
                if row["mrca_lineage"]:
                    info_dict[row["query"]]["scorpio_notes"] =  f'scorpio call: Alt alleles {row["alt_count"]}; Ref alleles {row["ref_count"]}; Amb alleles {row["ambig_count"]}; Oth alleles {row["other_count"]}'
                else:
                    info_dict[row["query"]]["scorpio_notes"] = ""
                    
        file_ending = fasta.split(".")[-1]
        if file_ending in ["gz","gzip","tgz"]:
            query = gzip.open(fasta, 'rt')
        else:
            query = open(fasta,"r")

        for l in query:
            if l[0]=='>':

                query_record = {}
                name = l.rstrip("\n").lstrip(">")
                modified_name=name.replace(" ","_").replace(",","_")
                if modified_name in name_dict:
                    query_record = info_dict[name_dict[modified_name]]
                    query_record["name"] = name
                    writer.writerow(query_record)
                else:
                    query_record = {
                        "name":name,
                        "hash":"",
                        "designated":"False",
                        "lineage":"Unassigned",
                        "qc_status":"fail",
                        "qc_notes":"failed to map",
                        "scorpio_mrca_lineage":"",
                        "scorpio_constellations":"",
                        "scorpio_incompatible_lineages":"",
                        "scorpio_support":"",
                        "scorpio_conflict":"",
                        "scorpio_notes":""
                    }








