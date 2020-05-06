from Bio import SeqIO
from Bio import Phylo
import sys
import os

if config.get("lineages_csv"):
    print("Going to run the global report summary")
else:
    config["lineages_csv"]=""


if config["lineages_csv"] != "":
    if config["write_tree"]==True:
        rule all:
            input:
                config["outfile"],
                config["outdir"] + "/global_lineage_information.csv",
                config["outdir"] + "/pangolin_trees/tree_file_names.txt"
    else: 
        rule all:
            input:
                config["outfile"],
                config["outdir"] + "/global_lineage_information.csv"
else:
    if config["write_tree"]==True:
        rule all:
            input:
                config["outfile"],
                config["outdir"] + "/pangolin_trees/tree_file_names.txt"
    else:
        rule all:
            input:
                config["outfile"]


rule expand_query_fasta:
    input:
        config["query_fasta"]
    params:
        config["query_sequences"]
    output:
        temp(expand(config["tempdir"] + '/{query}.fasta',query=config["query_sequences"]))
    run:
        for record in SeqIO.parse(input[0],"fasta"):
            with open(config["tempdir"] + f'/{record.id}.fasta',"w") as fw:
                fw.write(f">{record.id}\n{record.seq}\n")

rule profile_align_query:
    input:
        aln = config["representative_aln"],
        query = config["tempdir"] + '/{query}.fasta'
    output:
        temp(config["tempdir"] + "/{query}.aln.fasta")
    shell:
        "mafft --addprofile {input.query:q} {input.aln:q} > {output:q}"

rule iqtree_with_guide_tree:
    input:
        profile_aln = rules.profile_align_query.output,
        guide_tree = config["guide_tree"]
    output:
        temp(config["tempdir"] + "/{query}.aln.fasta.treefile"),
        temp(config["tempdir"] + "/{query}.aln.fasta.parstree"),
        temp(config["tempdir"] + "/{query}.aln.fasta.splits.nex"),
        temp(config["tempdir"] + "/{query}.aln.fasta.contree"),
        temp(config["tempdir"] + "/{query}.aln.fasta.log"),
        temp(config["tempdir"] + "/{query}.aln.fasta.ckp.gz"),
        temp(config["tempdir"] + "/{query}.aln.fasta.iqtree")
    run:
        iqtree_check = output[0].rstrip("treefile") + "iqtree"
        if os.path.exists(iqtree_check):
            print("Tree exists, going to rerun", iqtree_check)
            shell("iqtree -s {input.profile_aln:q} -bb 1000 -au -alrt 1000 -m HKY -g {input.guide_tree:q} -quiet -o 'outgroup_A' -redo")
        else:
            shell("iqtree -s {input.profile_aln:q} -bb 1000 -au -alrt 1000 -m HKY -g {input.guide_tree:q} -quiet -o 'outgroup_A'")

rule write_trees:
    input:
        trees = expand(config["tempdir"] + "/{query}.aln.fasta.treefile", query=config["query_sequences"]),
        key = config["key"]
    output:
        config["outdir"] + "/pangolin_trees/tree_file_names.txt"
    run:
        key_dict = {}
        with open(input.key, "r") as f:
            for l in f:
                l = l.rstrip('\n')
                taxon,key = l.split(",")
                key_dict[key] = taxon
        fout = open(output[0],"w")
        fout.write("taxon,treefile_name\n")
        for tree in input.trees:
            new_l = ""
            with open(tree, "r") as f:
                for l in f:
                    l = l.rstrip("\n")
                    new_l = l
                    for key in key_dict:
                        if key in l:
                            new_l = new_l.replace(key, key_dict[key])
                            taxon = key_dict[key].replace(" ","_").replace("/","_").replace("|","_")
                            fout.write(f"{key_dict[key]},{taxon}\n")
                            outdir = "/".join(output[0].split("/")[:-1])
                            with open(outdir + "/" + taxon + ".tree","w") as fw:
                                fw.write(new_l + "\n")

rule to_nexus:
    input:
        rules.iqtree_with_guide_tree.output
    output:
        temp(config["tempdir"] + "/{query}.nexus.tree")
    run:
        Phylo.convert(input[0], 'newick', output[0], 'nexus')

rule assign_lineage:
    input:
        tree = rules.to_nexus.output,
    params:
        query = "{query}",
        collapse=0.000005,

    output:
        temp(config["tempdir"] + "/{query}.txt")
    shell:
        """
        assign_lineage.py  --separator '_' --index 1 \
        --collapse_to_polytomies {params.collapse} --taxon '{params.query}' \
        --input {input.tree:q} --output {output:q}
        """
        
rule gather_reports:
    input:
        reports = expand(config["tempdir"] + "/{query}.txt", query=config["query_sequences"]),
        key=config["key"]
    params:
        version = config["lineages_version"]
    output:
        config["tempdir"] + "/lineage_report.pass_qc.csv"
    run:
        key_dict = {}
        with open(input.key, "r") as f:
            for l in f:
                l = l.rstrip('\n')
                taxon,key = l.split(",")
                key_dict[key] = taxon

        fw=open(output[0],"w")

        fw.write("taxon,lineage,SH-alrt,UFbootstrap,lineages_version,status,note\n")
        for lineage_report in input.reports:
            
            with open(lineage_report, "r") as f:
                for l in f:
                    l=l.rstrip()
                    key,lineage,alrt,bootstrap = l.split(",")
                    taxon = key_dict[key]
                    fw.write(f"{taxon},{lineage},{alrt},{bootstrap},{params.version},passed_qc,\n")
        fw.close()

rule add_failed_seqs:
    input:
        qcpass= config["tempdir"] + "/lineage_report.pass_qc.csv",
        qcfail= config["qc_fail"]
    params:
        version = config["lineages_version"]
    output:
        config["outfile"]
    run:
        fw = open(output[0],"w")
        with open(input.qcpass, "r") as f:
            for l in f:
                l=l.rstrip()
                fw.write(l + '\n')
        for record in SeqIO.parse(input.qcfail,"fasta"):
            desc_list = record.description.split(" ")
            note = ""
            for i in desc_list:
                if i.startswith("fail="):
                    note = i.lstrip("fail=")
            fw.write(f"{record.id},None,0,0,{params.version},fail,{note}\n")
        fw.close()

rule report_results:
    input:
        csv = config["outfile"],
        lineages_csv = config["lineages_csv"]
    output:
        config["outdir"] + "/global_lineage_information.csv"
    shell:
        """
        report_results.py \
        -p {input.csv} \
        -b {input.lineages_csv} \
        -o {output:q} 
        """