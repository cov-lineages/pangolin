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
        key = temp(config["tempdir"] + "/query_key.csv"),
        query_config = temp(config["tempdir"] + "/config.yaml")
    run:
        fkey = open(output.key, "w")
        ids = []
        c= 0
        
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                c+=1

                new_id = f"tax{params.pid}{c}tax"

                fkey.write(f"{record.id},{new_id}\n")
                fw.write(f">{new_id}\n{record.seq}\n")

                ids.append(new_id)
            
            print(f"{c} hashed sequences written")
        fkey.close()
        ids = ids.rstrip(',')
        with open(output.query_config, "w") as fconfig:
            fconfig.write("query_sequences: [")
            fconfig.write(ids + "]")
        query_sequence.store("query_store",c)


rule assign_lineages:
    input:
        snakefile = workflow.current_basedir+"/assign_query_lineage.smk",
        query = rules.pass_query_hash.output.fasta,
        key = rules.pass_query_hash.output.key,
        aln = rules.decrypt_aln.output,
        guide_tree = config["guide_tree"],
        query_config = config["tempdir"] + "/config.yaml"
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
        num_query_seqs = query_sequence.fetch("query_store")
        if num_query_seqs != 0:
            print(f"Passing {num_query_seqs} into processing pipeline.")
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "{params.force}"
                        "--directory {params.tempdir:q} "
                        "--configfile {input.query_config:q} "
                        "--config "
                        "outdir={params.outdir:q} "
                        "outfile={output.report:q} "
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


