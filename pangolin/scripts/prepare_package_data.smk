import csv
from Bio import SeqIO
import codecs
import collections

rule all:
    input:
        os.path.join(config["outdir"] , "anonymised.aln.fasta.treefile"),
        os.path.join(config["outdir"] , "anonymised.encrypted.aln.fasta"),
        os.path.join(config["outdir"] , "lineages.metadata.csv"),
        os.path.join(config["outdir"] , "defining_snps.csv"),
        os.path.join(config["outdir"] , "anonymised.aln.safe.fasta.treefile")

rule seqs_with_lineage:
    input:
        aln = config["fasta"],
        lineages = config["lineages"]
    output:
        fasta = os.path.join(config["outdir"], "sequences_with_lineage.fasta")
    run:
        to_write = []
        with open(input.lineages,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                to_write.append(row["taxon"])
        
        with open(output.fasta, "w") as fw:
            records = []
            for record in SeqIO.parse(input.aln,"fasta"):
                if "WH04" in record.id:
                    records.append(record)
                elif record.id in to_write:
                    records.append(record)
            SeqIO.write(records, fw, "fasta")

rule find_all_snps:
    input:
        aln = rules.seqs_with_lineage.output.fasta
    output:
        snps = os.path.join(config["outdir"] , "all_snps.csv")
    shell:
        """
        find_all_snps.py -a {input.aln:q} -o {output.snps:q}
        """

rule find_representatives:
    input:
        aln = rules.seqs_with_lineage.output.fasta,
        lineages = config["lineages"],
        snps = rules.find_all_snps.output.snps,
        metadata = config["metadata"]
    output:
        reps = os.path.join(config["outdir"] , "representative_seqs.csv"),
        defining = os.path.join(config["outdir"] , "defining_snps.csv"),
        mask = os.path.join(config["outdir"] , "to_mask.csv")
    shell:
        """categorise_snps.py -a {input.aln:q} -l {input.lineages:q} \
                --snps {input.snps:q} \
                --metadata {input.metadata:q} \
                --representative-seqs-out {output.reps:q} \
                --defining-snps-out {output.defining:q} \
                --mask-out {output.mask:q} 
            """

rule extract_representative_sequences:
    input:
        aln = rules.seqs_with_lineage.output.fasta,
        lineages = config["lineages"],
        metadata = config["metadata"],
        mask = os.path.join(config["outdir"] , "to_mask.csv"),
        representatives = os.path.join(config["outdir"] , "representative_seqs.csv")
    output:
        representatives = os.path.join(config["outdir"] , "representative_sequences.fasta"),
        metadata = os.path.join(config["outdir"] , "lineages.metadata.csv")
    shell:
        """
        get_masked_representatives.py \
            -a {input.aln} \
            -l {input.lineages} \
            --representatives {input.representatives} \
            --to-mask {input.mask} \
            --metadata {input.metadata} \
            --representative-seqs-out {output.representatives} \
            --metadata-out {output.metadata} 
        """

rule mafft_representative_sequences:
    input:
        rules.extract_representative_sequences.output.representatives
    threads: workflow.cores
    params:
        cores = workflow.cores
    output:
        os.path.join(config["outdir"] , "representative_sequences.aln.fasta")
    shell:
        "mafft --thread {params.cores} {input[0]:q} > {output[0]:q}"

rule anonymise_headers:
    input:
        rules.mafft_representative_sequences.output
    output:
        fasta = os.path.join(config["outdir"] , "anonymised.aln.fasta"),
        key = os.path.join(config["outdir"] , "tax_key.csv")
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
    threads: workflow.cores
    params:
        cores = workflow.cores
    output:
        os.path.join(config["outdir"] , "anonymised.aln.fasta.treefile")
    shell:
        "iqtree -s {input[0]:q} -nt AUTO -bb 10000 -m HKY -redo -au -alrt 1000 -o 'outgroup_A'"


rule anonymise_headers_safe:
    input:
        rules.mafft_representative_sequences.output
    output:
        fasta = os.path.join(config["outdir"] , "anonymised.aln.safe.fasta"),
        key = os.path.join(config["outdir"] , "tax_key.safe.csv")
    run:
        fkey = open(output.key, "w")
        fkey.write("taxon,key\n")
        key_dict = {}
        c = 0
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                c+=1
                name,lineage = record.id.split('|')
                lineage_list = lineage.split(".")
                if lineage_list[-1].startswith("p"):
                    lineage = ".".join(lineage_list[:-1])
                new_id = ""
                if "WH04" in record.id: # this is the WH04 GISAID ID
                    print(record.id)
                    new_id = f"outgroup_A"
                else:
                    new_id = f"{c}_{lineage}"

                fkey.write(f"{record.id},{new_id}\n")
                fw.write(f">{new_id}\n{record.seq}\n")
            
            print(f"{c+1} safe anonymised sequences written to {output.fasta}")
        fkey.close()

rule encrypt_fasta_safe:
    input:
        rules.anonymise_headers_safe.output.fasta
    output:
        os.path.join(config["outdir"] , "anonymised.encrypted.aln.safe.fasta")
    run:
        c = 0
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                encrypted_seq = codecs.encode(str(record.seq), 'rot_13')
                fw.write(f">{record.id}\n{encrypted_seq}\n")
                c+=1
        print(f"Encrypted {c} sequences with rot13")

rule iqtree_representative_sequences_safe:
    input:
        rules.anonymise_headers_safe.output.fasta
    threads: workflow.cores
    params:
        cores = workflow.cores
    output:
        os.path.join(config["outdir"] , "anonymised.aln.fasta.safe.treefile")
    shell:
        "iqtree -s {input[0]:q} -nt AUTO -bb 10000 -m HKY -redo -au -alrt 1000 -o 'outgroup_A'"

