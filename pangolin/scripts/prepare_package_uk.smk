import csv
from Bio import SeqIO
import codecs

# requires metadata and fasta file in config

rule all:
    input:
        os.path.join(config["outdir"] , "anonymised.aln.fasta.treefile"),
        os.path.join(config["outdir"] , "anonymised.encrypted.aln.fasta")

rule extract_representative_sequences:
    input:
        metadata = config["metadata"],
        fasta = config["fasta"]
    output:
        os.path.join(config["outdir"] , "representative_sequences.fasta")
    run:
        tax_dict = {}
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["representative"] == '1':
                    tax_dict[row["name"]] = row["UK_cluster"]
        print(f"Number of seqs in tax dict is {len(tax_dict)}")
        fw = open(output[0], "w")
        c = 0
        for record in SeqIO.parse(input.fasta,"fasta"):
            if record.id in tax_dict:
                c+=1
                record_list = record.id
                new_id = f"{record.id}|{tax_dict[record.id]}"
                fw.write(f">{new_id}\n{record.seq}\n")
        print(f"{c} representative sequences written to {output[0]}")
        fw.close()

rule mafft_representative_sequences:
    input:
        rules.extract_representative_sequences.output
    output:
        os.path.join(config["outdir"] , "representative_sequences.aln.fasta")
    shell:
        "mafft {input[0]:q} > {output[0]:q}"

rule anonymise_headers:
    input:
        rules.mafft_representative_sequences.output
    output:
        fasta = os.path.join(config["outdir"] , "anonymised.aln.fasta"),
        key = os.path.join(config["outdir"], "tax_key.csv")
    run:
        fkey = open(output.key, "w")
        fkey.write("taxon,key\n")
        key_dict = {}
        c = 0
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                c+=1
                record_list = record.id.split('|')
                lineage = record_list[1]
                new_id = ""
                if "WH04" in record.id: # this is the WH04 GISAID ID
                    print(record.id)
                    new_id = f"outgroup_NonUK"
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
        os.path.join(config["outdir"] , "anonymised.encrypted.aln.fasta")
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
        os.path.join(config["outdir"] , "anonymised.aln.fasta.treefile")
    shell:
        "iqtree -s {input[0]:q} -bb 1000 -m HKY -o 'outgroup_NonUK' -redo"
