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
        os.path.join(config["outdir"] , "anonymised.aln.safe.fasta.treefile"),
        os.path.join(config["outdir"] , "anonymised.encrypted.aln.safe.fasta")

rule phylotype_to_metadata:
    input:
        phylotype = config["phylotypes"],
        metadata = config["metadata"]
    output:
        os.path.join(config["outdir"], "metadata_with_phylotypes.csv")
    run:
        phylotype_dict = {}
        with open(input.phylotype,"r") as f:
            for l in f:
                l = l.rstrip("\n")
                name,phylotype=l.split(',')
                phylotype_dict[name]=phylotype
            
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)

            header = reader.fieldnames
            header.append("phylotype")
            with open(output[0], "w") as fw:
                writer = csv.DictWriter(fw, fieldnames=header,lineterminator='\n')
                writer.writeheader()
                for row in reader:
                    if row["sequence_name"] in phylotype_dict:
                        new_row = row
                        name = new_row["sequence_name"]
                        phylotype = phylotype_dict[name]
                        new_row["phylotype"] = phylotype

                        writer.writerow(new_row)

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
                if "Wuhan/WH04/2020" == record.id:
                    records.append(record)
                elif record.id in to_write:
                    records.append(record)
            SeqIO.write(records, fw, "fasta")

#input config all snp file

rule find_all_snps:
    input:
        aln = rules.seqs_with_lineage.output.fasta
    output:
        snps = os.path.join(config["outdir"] , "all_snps.csv")
    shell:
        """
        find_all_snps.py -a {input.aln:q} -o {output.snps:q}
        """

rule find_polytomies:
    input:
        tree = config["global_tree"]
    output:
        outfile = os.path.join(config["outdir"] , "all_polytomies.csv")
    shell:
        """
        get_polytomy.py \
            --global-tree {input.tree} \
            --outfile {output.outfile}
        """

rule find_basal_polytomies:
    input:
        polytomies = rules.find_polytomies.output.outfile,
        lineages = config["lineages"],
        metadata = os.path.join(config["outdir"], "metadata_with_phylotypes.csv")
    output:
        outfile = os.path.join(config["outdir"] , "basal_polytomy_taxa.csv")
    shell:
        """
        get_basal_polytomy.py \
            --polytomies {input.polytomies} \
            --lineages {input.lineages} \
            --metadata {input.metadata} \
            --outfile {output.outfile}
        """

rule find_representatives:
    input:
        aln = rules.seqs_with_lineage.output.fasta,
        lineages = config["lineages"],
        tree = config["global_tree"],
        snps = os.path.join(config["outdir"] , "all_snps.csv"), #config["all_snps"],
        include = config["to_include"],
        polytomies = rules.find_basal_polytomies.output.outfile
    output:
        reps = os.path.join(config["outdir"] , "representative_seqs.csv"),
        defining = os.path.join(config["outdir"] , "defining_snps.csv"),
        mask = os.path.join(config["outdir"] , "singletons.csv"),
    shell:
        """categorise_snps.py \
                -a {input.aln:q} \
                -l {input.lineages:q} \
                --snps {input.snps:q} \
                --global-tree {input.tree} \
                --polytomy {input.polytomies:q} \
                --to-include {input.include} \
                --representative-seqs-out {output.reps:q} \
                --defining-snps-out {output.defining:q} \
                --mask-out {output.mask:q} \
            """

rule extract_representative_sequences:
    input:
        aln = rules.seqs_with_lineage.output.fasta,
        lineages = config["lineages"],
        metadata = os.path.join(config["outdir"], "metadata_with_phylotypes.csv"),
        mask = rules.find_representatives.output.mask,
        representatives = rules.find_representatives.output.reps
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
                if "Wuhan/WH04/2020" in record.id: # this is the WH04 GISAID ID
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
        cores = workflow.cores/2
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
                if "Wuhan/WH04/2020" in record.id: # this is the WH04 GISAID ID
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
        cores = workflow.cores/2
    output:
        os.path.join(config["outdir"] , "anonymised.aln.safe.fasta.treefile")
    shell:
        "iqtree -s {input[0]:q} -nt AUTO -bb 10000 -m HKY -redo -au -alrt 1000 -o 'outgroup_A'"

