import csv
from Bio import SeqIO
import codecs
import collections

rule all:
    input:
        config["outdir"] + "/anonymised.aln.fasta.treefile",
        config["outdir"] + "/anonymised.encrypted.aln.fasta",
        config["outdir"] + "/lineages.metadata.csv",
        defining = config["outdir"] + "/defining_snps.csv"

rule find_representatives:
    input:
        aln = config["fasta"],
        lineages = config["lineages"]
    output:
        snps = config["outdir"] + "/representative_seqs.csv",
        defining = config["outdir"] + "/defining_snps.csv",
        mask = config["outdir"] + "/to_mask.csv"
    shell:
        """all_snps.py -a {input.aln:q} -l {input.lineages:q} \
                --all_snps {output.snps:q} \
                --defining_snps {output.defining:q} \
                --to_mask {output.mask:q} 
            """

rule extract_representative_sequences:
    input:
        aln = config["fasta"],
        mask = config["outdir"] + "/to_mask.csv",
        metadata = config["metadata"],
        lineages = config["lineages"],
        representatives = config["outdir"] + "/representative_seqs.csv"
    output:
        representatives = config["outdir"] + "/representative_sequences.fasta",
        metadata = config["outdir"] + "/lineages.metadata.csv",
    shell:
        """
        get_masked_representatives.py \
            -r {input.representatives} \
            -m {input.mask} \
            -l {input.lineages} \
            -a {input.aln} \
            --metadata {input.metadata} \
            -o {output.representatives} \
            --metadata-out {output.metadata} 
        """

rule mafft_representative_sequences:
    input:
        rules.extract_representative_sequences.output.representatives
    output:
        config["outdir"] + "/representative_sequences.aln.fasta"
    shell:
        "mafft {input[0]:q} > {output[0]:q}"

rule anonymise_headers:
    input:
        rules.mafft_representative_sequences.output
    output:
        fasta = config["outdir"] + "/anonymised.aln.fasta",
        key = config["outdir"] + "/tax_key.csv"
    run:
        fkey = open(output.key, "w")
        fkey.write("taxon,key\n")
        key_dict = {}
        c = 0
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                c+=1
                name,lineage = record.id.split('|')
                new_id = ""
                if "WH04" in record.id: # this is the WH04 GISAID ID
                    print(record.id)
                    new_id = f"outgroup_A"
                else:
                    new_id = f"{c}_{lineage}"

                fkey.write(f"{record.id},{new_id}\n")
                fw.write(f">{new_id}\n{record.seq}\n")
            
            print(f"{c+1} anonymised sequences written to {output.fasta}")
        fkey.close()

rule encrypt_fasta:
    input:
        rules.anonymise_headers.output.fasta
    output:
        config["outdir"] + "/anonymised.encrypted.aln.fasta"
    run:
        c = 0
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                encrypted_seq = codecs.encode(str(record.seq), 'rot_13')
                fw.write(f">{record.id}\n{encrypted_seq}\n")
                c+=1
        print(f"Encrypted {c} sequences with rot13")

rule iqtree_representative_sequences:
    input:
        rules.anonymise_headers.output.fasta
    output:
        config["outdir"] + "/anonymised.aln.fasta.treefile"
    shell:
        "iqtree -s {input[0]:q} -bb 10000 -m HKY -redo -au -alrt 1000 -o 'outgroup_A'"
