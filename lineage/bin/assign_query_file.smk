query_sequence = PersistentDict("query_store")


rule iqtree_safe_query:
    input:
        config["query_fasta"]
    output:
        temp(config["outdir"] + "/temp/query_fasta.iqtree_friendly.fasta")
    run:
        fw = open(output[0],"w")
        for record in SeqIO.parse(input[0],"fasta"):
            new_id = record.id.replace("/","___")
            fw.write(f">{new_id}\n{record.seq}\n")
        fw.close()

rule assess_queries:
    input:
        rules.iqtree_safe_query.output
    output:
        t = temp(config["outdir"] + "/temp/temp_report.txt")
    run:
        fw = open(output[0],"w")
        ids = ''
        for record in SeqIO.parse(input[0],"fasta"):
            gisaid_id = record.id.split("|")[1]
            if gisaid_id != "":
                ids += gisaid_id + ","
                
                fw.write(f"{gisaid_id},{record.id}\n")
            else:
                cog_id = record.id.split("|")[0].split('___')[1]
                ids += cog_id + ","
        fw.close()
        ids = ids.rstrip(',')
        query_sequence.store("query_store",ids)

rule process_sample:
    input:
        rules.assess_queries.output.t,
        config=workflow.current_basedir+"/../config.yaml",
        snakefile = workflow.current_basedir+"/assign_query_lineage.smk",
        query = rules.iqtree_safe_query.output,
        aln = config["representative_aln"],
        guide_tree = config["guide_tree"]
    params:
        outdir= config["outdir"],
        path = workflow.current_basedir,
    output:
        report = config["outdir"] + "/lineage_report.txt"
    run:
        query_sequences = query_sequence.fetch("query_store")
        if query_sequences != "":
            print("Passing {} into processing pipeline.".format(query_sequences))
            config["query_sequences"]= query_sequences
            shell("snakemake --nolock --snakefile {input.snakefile:q} "
                        "--configfile {input.config:q} "
                        "--config "
                        "query_sequences={config[query_sequences]} "
                        "outdir={params.outdir:q} "
                        "query_fasta={input.query:q} "
                        "representative_aln={input.aln:q} "
                        "guide_tree={input.guide_tree:q} "
                        "--rerun-incomplete --cores 1")
        else:
            shell("touch {output.report}")



