query_sequence = PersistentDict("query_store")

rule decrypt_aln:
    input:
        config["representative_aln"]
    output:
        config["outdir"] + "/temp/anonymised.aln.fasta"
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
    output:
        t = temp(config["outdir"] + "/temp/temp.txt"),
        fasta = config["outdir"] + "/temp/query.fasta",
        key = config["outdir"] + "/temp/query_key.csv"
    run:
        fkey = open(output.key, "w")
        ids = ''
        c= 0
        with open(output.fasta, "w") as fw:
            for record in SeqIO.parse(input[0],"fasta"):
                c+=1
                record_list = record.id.split('|')
                lineage = record_list[2]
                
                new_id = f"{c}_{lineage}"

                fkey.write(f"{record.id},{new_id}\n")
                fw.write(f">{new_id}\n{record.seq}\n")

                ids+=new_id + ','
            
            print(f"{c+1} hashed sequences written")
        fkey.close()
        
        ids = ids.rstrip(',')
        query_sequence.store("query_store",ids)
        shell("touch {output.t}")


rule assign_lineages:
    input:
        rules.pass_query_hash.output.t,
        config=workflow.current_basedir+"/../config.yaml",
        snakefile = workflow.current_basedir+"/assign_query_lineage.smk",
        query = rules.pass_query_hash.output.fasta,
        key = rules.pass_query_hash.output.key,
        aln = rules.decrypt_aln.output,
        guide_tree = config["guide_tree"]
    params:
        outdir= config["outdir"],
        path = workflow.current_basedir,
        cores = workflow.cores
    output:
        report = config["outdir"] + "/lineage_report.csv"
    run:
        query_sequences = query_sequence.fetch("query_store")
        num_query_seqs = len(query_sequences.split(","))
        if query_sequences != "":
            print(f"Passing {num_query_seqs} into processing pipeline.")
            config["query_sequences"]= query_sequences
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "--configfile {input.config:q} "
                        "--config "
                        "query_sequences={config[query_sequences]} "
                        "outdir={params.outdir:q} "
                        "query_fasta={input.query:q} "
                        "representative_aln={input.aln:q} "
                        "guide_tree={input.guide_tree:q} "
                        "key={input.key} "
                        "--cores {params.cores}")
        else:
            shell("touch {output.report}")



