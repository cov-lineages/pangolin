from Bio import SeqIO
from Bio import Phylo

from pytools.persistent_dict import PersistentDict

if config.get("outdir"):
    config["outdir"] = config["outdir"].rstrip("/")
else:
    config["outdir"] = "analysis"


taxon_dict = PersistentDict("taxon")

config["query_sequences"]=[i for i in config["query_sequences"].split(',')]

rule all:
    input:
        expand(config["outdir"] + "/temp/query_alignments/{query}.aln.fasta.treefile", query=config["query_sequences"]),
        config["outdir"] + "/lineage_report.txt"

rule expand_query_fasta:
    input:
        config["query_fasta"]
    params:
        config["query_sequences"]
    output:
        config["outdir"] + '/temp/temp.txt')
        expand(config["outdir"] + '/temp/expanded_query/{query}.fasta',query=config["query_sequences"])
    run:
        shell("touch {output[0]:q}")
        for record in SeqIO.parse(input[0],"fasta"):
            gisaid_id = record.id.split("|")[1]
            print(record.id, gisaid_id)
            if gisaid_id != "":
                taxon_dict.store(gisaid_id, record.id)
                print(gisaid_id)
                with open(config["outdir"] + f'/temp/expanded_query/{gisaid_id}.fasta',"w") as fw:
                    fw.write(f">{record.id}\n{record.seq}\n")
            else:
                cog_id = record.id.split("|")[0].split('___')[1]
                print(cog_id)
                taxon_dict.store(cog_id, record.id)
                with open(config["outdir"] + f'/temp/expanded_query/{cog_id}.fasta',"w") as fw:
                    fw.write(f">{record.id}\n{record.seq}\n")
 
rule profile_align_new_fasta:
    input:
        aln = config["representative_aln"],
        query = config["outdir"] + '/temp/expanded_query/{query}.fasta'
    output:
        config["outdir"] + "/temp/query_alignments/{query}.aln.fasta"
    shell:
        "mafft-profile {input.aln:q} {input.query:q} > {output:q}"

rule iqtree_with_guide_tree:
    input:
        profile_aln = rules.profile_align_new_fasta.output,
        guide_tree = config["guide_tree"]
    output:
        config["outdir"] + "/temp/query_alignments/{query}.aln.fasta.treefile"
    shell:
        "iqtree -s {input.profile_aln:q} -bb 1000 -m HKY -g {input.guide_tree:q}"

rule iqtree_to_nexus:
    input:
        rules.iqtree_with_guide_tree.output
    output:
        config["outdir"] + "/temp/query_alignments/{query}.aln.fasta.nexus.tree"
    run:
        Phylo.convert(input[0], 'newick', output[0], 'nexus')

rule assign_lineage:
    input:
        tree = rules.iqtree_to_nexus.output,
    params:
        query = "{query}"
    output:
        config["outdir"] + "/temp/reports/{query}.csv"
    run:
        taxon = taxon_dict.fetch(params.query)
        print(params.query, taxon)
        shell_start = f"clusterFunk subtype  --separator '|' --index 2 --collapse_to_polytomies --taxon '{taxon}'"
        shell(shell_start + " --input {input.tree:q} --output {output:q}")
        
rule gather_reports:
    input:
        expand(config["outdir"] + "/temp/reports/{query}.csv", query=config["query_sequences"])
    output:
        config["outdir"] + "/lineage_report.txt"
    run:
        fw=open(output[0],"w")
        fw.write("taxon,lineage\n")
        for lineage_report in input:
            with open(lineage_report, "r") as f:
                for l in f:
                    l=l.rstrip()
                    l = l.replace("___","/")
                    fw.write(l+'\n')
        fw.close()

