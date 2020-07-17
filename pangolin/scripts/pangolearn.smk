#!/usr/bin/env python

# config["trained_model"] = trained_model
import csv
from Bio import SeqIO
import os

##### Configuration #####

if config.get("trained_model"):
    config["trained_model"] = os.path.join(workflow.current_basedir,'..', config["trained_model"])

if config.get("header_file"):
    config["header_file"] = os.path.join(workflow.current_basedir,'..', config["header_file"])

if config.get("lineages_csv"):
    lineages_csv_path = os.path.join(workflow.current_basedir,'..', config["lineages_csv"])
    config["lineages_csv"]=f"lineages_csv={lineages_csv_path} "
else:
    config["lineages_csv"]=""

if config.get("force"):
    config["force"] = "--forceall "
else:
    config["force"] = ""

if config.get("lineages_csv"):
    print("Going to run the global report summary")
else:
    config["lineages_csv"]=""

config["pid"] = str(os.getpid())

##### Target rules #####

if config["lineages_csv"] != "":
    rule all:
        input:
            config["outfile"],
            os.path.join(config["outdir"],"global_lineage_information.csv")
else:
    rule all:
        input:
            config["outfile"]

rule minimap2_to_reference:
    input:
        fasta = config["query_fasta"],
        reference = config["reference_fasta"]
    output:
        sam = os.path.join(config["tempdir"],"reference_mapped.sam")
    shell:
        """
        minimap2 -a -x asm5 {input.reference:q} {input.fasta:q} > {output.sam:q}
        """

rule datafunk_trim_and_pad:
    input:
        sam = rules.minimap2_to_reference.output.sam,
        reference = config["reference_fasta"]
    params:
        trim_start = config["trim_start"],
        trim_end = config["trim_end"],
        insertions = os.path.join(config["tempdir"],"insertions.txt")
    output:
        fasta = os.path.join(config["tempdir"],"post_qc_query.aligned.fasta")
    shell:
        """
        datafunk sam_2_fasta \
          -s {input.sam:q} \
          -r {input.reference:q} \
          -o {output.fasta:q} \
          -t [{params.trim_start}:{params.trim_end}] \
          --pad \
          --log-inserts 
        """

rule pangolearn:
    input:
        fasta = rules.datafunk_trim_and_pad.output.fasta,
        model = config["trained_model"],
        header = config["header_file"]
    output:
        os.path.join(config["tempdir"],"lineage_report.pass_qc.csv")
    shell:
        # should output a csv file with no headers but with columns similar to:
        # "taxon,lineage,SH-alrt,UFbootstrap"
        """
        pangolearn.py --header-file {input.header} --model-file {input.model} --fasta {input.fasta} -o {output[0]}
        """

rule add_failed_seqs:
    input:
        qcpass= os.path.join(config["tempdir"],"lineage_report.pass_qc.csv"),
        qcfail= config["qc_fail"]
    params:
        version = config["lineages_version"]
    output:
        config["outfile"]
    run:
        fw = open(output[0],"w")
        fw.write("taxon,lineage,SH-alrt,UFbootstrap,lineages_version,status,note\n")

        with open(input.qcpass, "r") as f:
            for l in f:
                l=l.rstrip('\n')
                fw.write(f"{l},{params.version},passed_qc,\n")

        for record in SeqIO.parse(input.qcfail,"fasta"):
            desc_list = record.description.split(" ")
            note = ""
            for i in desc_list:
                if i.startswith("fail="):
                    note = i.lstrip("fail=")
            # needs to mirror the structure of the output from pangolearn
            fw.write(f"{record.id},None,0,0,{params.version},fail,{note}\n")

        fw.close()

rule report_results:
    input:
        csv = config["outfile"],
        lineages_csv = config["lineages_csv"]
    output:
        os.path.join(config["outdir"],"global_lineage_information.csv")
    shell:
        """
        report_results.py \
        -p {input.csv} \
        -b {input.lineages_csv} \
        -o {output:q} 
        """
