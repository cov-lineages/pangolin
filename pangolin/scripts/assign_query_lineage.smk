from Bio import SeqIO
from Bio import Phylo
import sys

config["query_sequences"]=[i for i in config["query_sequences"].split(',')]

if config.get("lineages_csv"):
    rule all:
        input:
            config["outdir"] + "/lineage_report.csv",
            config["outdir"] + "/global_lineage_information.csv"
else:
    rule all:
        input:
            config["outdir"] + "/lineage_report.csv"


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
        "mafft-profile {input.aln:q} {input.query:q} > {output:q}"

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
        query = "{query}"
    output:
        temp(config["tempdir"] + "/{query}.txt")
    shell:
        """
        assign_lineage.py  --separator '_' --index 1 \
        --collapse_to_polytomies --taxon '{params.query}' \
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
                    tokens = l.split(",")
                    lineage,support = tokens[1],tokens[2]
                    taxon = key_dict[tokens[0]]
                    bootstrap = ""
                    print(support)
                    support = support.split("/")
                    if len(support) == 4: 
                        old_alrt,old_bs,alrt,ufboot = support
                        bootstrap = ufboot.split('.')[0]
                        alrt = alrt.split('.')[0]
                        print("alrt",alrt,"bootstrap",bootstrap)
                    elif len(support) == 2:
                        alrt,ufboot = support
                        bootstrap = ufboot.split('.')[0]
                        alrt = alrt.split('.')[0]
                        print("alrt",alrt,"bootstrap",bootstrap)
                    else:
                        alrt=0
                        bootstrap=0

                    fw.write(f"{taxon},{lineage},{alrt},{bootstrap},{params.version},passed_qc,\n")
        fw.close()

rule add_failed_seqs:
    input:
        qcpass= config["tempdir"] + "/lineage_report.pass_qc.csv",
        qcfail= config["qc_fail"]
    params:
        version = config["lineages_version"]
    output:
        config["outdir"] + "/lineage_report.csv"
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
        csv = config["outdir"] + "/lineage_report.csv",
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