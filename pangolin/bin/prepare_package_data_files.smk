
rule extract_representative_sequences:
    input:
        metadata = config["lineage_metadata"],
        fasta = config["fasta"]
    output:
        config["outdir"] + "/temp/representative_sequences.fasta"
    run:
        tax_dict = {}
        with open(input.metadata,newline="") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["representative"] == '1':
                    tax_dict[row["header"]] = row["lineage"]
        fw = open(output[0], "w")
        c = 0
        for record in SeqIO.parse(input.fasta,"fasta"):
            if record.id in tax_dict:
                c+=1
                
                record_list = record.id.split('|')
                country = record_list[0].split("/")[0]
                record_list[2] = tax_dict[record.id]
                new_id = '|'.join(record_list)
                fw.write(f">{new_id}\n{record.seq}\n")
        print(f"{c} representative sequences written to {output[0]}")
        fw.close()

rule mafft_representative_sequences:
    input:
        config["outdir"] + "/temp/representative_sequences.fasta"
    output:
        config["outdir"] + "/temp/representative_sequences.aln.fasta"
    shell:
        "mafft {input[0]:q} > {output[0]:q}"

rule iqtree_safe_headers:
    input:
        config["outdir"] + "/temp/representative_sequences.aln.fasta"
    output:
        config["outdir"] + "/temp/representative_sequences.aln.iqtree_friendly.fasta"
    run:
        fw = open(output[0],"w")
        for record in SeqIO.parse(input[0],"fasta"):
            new_id = record.id.replace("/","___")
            fw.write(f">{new_id}\n{record.seq}\n")
        fw.close()
            

rule iqtree_representative_sequences:
    input:
        config["outdir"] + "/temp/representative_sequences.aln.iqtree_friendly.fasta"
    output:
        config["outdir"] + "/temp/representative_sequences.aln.iqtree_friendly.fasta.treefile"
    shell:
        "iqtree -s {input[0]:q} -bb 1000 -m HKY -o 'Wuhan___WH04___2020|EPI_ISL_406801|A|Wuhan|||2020-01-05'"

rule iqtree_restore_headers:
    input:
        config["outdir"] + "/temp/representative_sequences.aln.iqtree_friendly.fasta.treefile"
    output:
        config["outdir"] + "/temp/representative_sequences.aln.fasta.treefile"
    run:
        fw = open(output[0],"w")
        with open(input[0],"r") as f:
            for l in f:
                l = l.rstrip("\n")
                l = l.replace("___","/")
            fw.write(f"{l}\n")
        fw.close()
