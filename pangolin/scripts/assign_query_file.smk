query_sequence = PersistentDict("query_store")

rule decrypt_aln:
    input:
        config["representative_aln"]
    output:
        temp(config["tempdir"] +"/anonymised.aln.fasta")
    run:
        c = 0
        with open(output[0],"w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                decrypted_seq = codecs.decode(str(record.seq), 'rot_13')
                fw.write(f">{record.id}\n{decrypted_seq}\n")
                c+=1
        print(f"Decrypted {c} sequences")

rule pass_query_hash:
    input:
        config["query_fasta"]
    params:
        pid = config["pid"]
    output:
        fasta = temp(config["tempdir"] + "/query.fasta"),
        key = temp(config["tempdir"] + "/query_key.csv")
    run:
        fkey = open(output.key, "w")
        ids = ''
        c= 0
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                c+=1

                new_id = f"tax{params.pid}{c}tax"

                fkey.write(f"{record.id},{new_id}\n")
                fw.write(f">{new_id}\n{record.seq}\n")

                ids+=new_id + ','
            
            print(f"{c+1} hashed sequences written")
        fkey.close()
        
        ids = ids.rstrip(',')
        query_sequence.store("query_store",ids)


rule assign_lineages:
    input:
        snakefile = workflow.current_basedir+"/assign_query_lineage.smk",
        query = rules.pass_query_hash.output.fasta,
        key = rules.pass_query_hash.output.key,
        aln = rules.decrypt_aln.output,
        guide_tree = config["guide_tree"]
    params:
        outdir= config["outdir"],
        tempdir= config["tempdir"],
        qcfail=config["qc_fail"],
        path = workflow.current_basedir,
        cores = workflow.cores,
        force = config["force"],
        write_tree=config["write_tree"],
        lineages_csv=config["lineages_csv"],
        version=config["lineages_version"]
    output:
        report = config["outfile"]
    run:
        query_sequences = query_sequence.fetch("query_store")
        num_query_seqs = len(query_sequences.split(","))
        if query_sequences != "":
            print(f"Passing {num_query_seqs} into processing pipeline.")
            config["query_sequences"]= query_sequences
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force}"
                        "--directory {params.tempdir:q} "
                        "--config "
                        "query_sequences={config[query_sequences]} "
                        "outdir={params.outdir:q} "
                        "outfile={output.report} "
                        "write_tree={params.write_tree} "
                        "tempdir={params.tempdir:q} "
                        "query_fasta={input.query:q} "
                        "qc_fail={params.qcfail:q} "
                        "representative_aln={input.aln:q} "
                        "lineages_version={params.version} "
                        "{params.lineages_csv}"
                        "guide_tree={input.guide_tree:q} "
                        "key={input.key:q} "
                        "--cores {params.cores}")
        else:
            fw = open(output.report,"w")
            fw.write("taxon,lineage,SH-alrt,UFbootstrap,lineages_version,status,note\n")
            for record in SeqIO.parse(params.qcfail,"fasta"):
                desc_list = record.description.split(" ")
                note = ""
                for i in desc_list:
                    if i.startswith("fail="):
                        note = i.lstrip("fail=")
                fw.write(f"{record.id},None,0,0,{params.version},fail,{note}\n")
            fw.close()


