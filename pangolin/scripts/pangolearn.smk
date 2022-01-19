#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
import gzip
from pangolin.utils.log_colours import green,cyan,red
from pangolin.utils.hash_functions import get_hash_string
from pangolin.utils.sequence_qc import sequence_qc

import pangolin.pangolearn.pangolearn as pangolearn

##### Configuration #####

if config.get("trained_model"):
    config["trained_model"] = os.path.join(workflow.current_basedir,'..', config["trained_model"])

if config.get("header_file"):
    config["header_file"] = os.path.join(workflow.current_basedir,'..', config["header_file"])

##### Utility functions #####


##### Report options #####
UNASSIGNED_LINEAGE_REPORTED="None"

##### Target rules #####


rule all:
    input:
        config["outfile"],
        os.path.join(config["tempdir"],"VOC_report.scorpio.csv")


rule align_to_reference:
    input:
        fasta = config["query_fasta"],
        reference = config["reference_fasta"]
    params:
        trim_start = 265,
        trim_end = 29674,
        sam = os.path.join(config["tempdir"],"mapped.sam")
    output:
        fasta = os.path.join(config["aligndir"],"sequences.aln.fasta")
    log:
        os.path.join(config["tempdir"], "logs/minimap2_sam.log")
    shell:
        """
        minimap2 -a -x asm20 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} '{input.fasta}' -o {params.sam:q} &> {log:q} 
        gofasta sam toMultiAlign \
            -s {params.sam:q} \
            -t {workflow.cores} \
            --reference {input.reference:q} \
            --trimstart {params.trim_start} \
            --trimend {params.trim_end} \
            --trim \
            --pad > '{output.fasta}'
        """


rule hash_sequence_assign:
    input:
        fasta = rules.align_to_reference.output.fasta
    output:
        designated = os.path.join(config["tempdir"],"hash_assigned.csv"),
        for_inference = os.path.join(config["tempdir"],"not_assigned.fasta")
    params:
        skip_designation_hash = config["skip_designation_hash"]
    run:
        set_hash = {}
        with open(config["designated_hash"],"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                set_hash[row["seq_hash"]] = row["lineage"]
        
        with open(output.designated,"w") as fw:
            fw.write("taxon,lineage\n")
            with open(output.for_inference, "w") as fseq:
                for record in SeqIO.parse(input.fasta, "fasta"):
                    if record.id!="reference":
                        hash_string = get_hash_string(record)
                        if not params.skip_designation_hash and hash_string in set_hash:
                            fw.write(f"{record.id},{set_hash[hash_string]}\n")
                        else:
                            fseq.write(f">{record.description}\n{record.seq}\n")


rule scorpio:
    input:
        fasta = rules.align_to_reference.output.fasta,
    params:
        constellation_files = " ".join(config["constellation_files"])
    output:
        report = os.path.join(config["tempdir"],"VOC_report.scorpio.csv")
    threads:
        workflow.cores
    log:
        os.path.join(config["tempdir"], "logs/scorpio.log")
    shell:
        """
        scorpio classify \
        -i {input.fasta:q} \
        -o {output.report:q} \
        -t {workflow.cores} \
        --output-counts \
        --constellations {params.constellation_files} \
        --pangolin \
        --list-incompatible \
        --long &> {log:q}
        """

rule get_constellations:
    params:
        constellation_files = " ".join(config["constellation_files"])
    output:
        list = os.path.join(config["tempdir"], "get_constellations.txt")
    shell:
        """
        scorpio list \
        --constellations {params.constellation_files} \
        --pangolin > {output.list:q}
        """


rule sequence_qc:
    input:
    output:
    run:


rule pangolearn:
    input:
        fasta = rules.cache_sequence_assign.output.for_inference,
        model = config["trained_model"],
        header = config["header_file"],
        reference = config["reference_fasta"]
    output:
        os.path.join(config["tempdir"],"lineage_report.pass_qc.csv")
    run:
        pangolearn.assign_lineage(input.header,input.model,input.reference,input.fasta,output[0])

rule add_failed_seqs:
    input:
        qcpass= os.path.join(config["tempdir"],"lineage_report.pass_qc.csv"),
        qcfail= config["qc_fail"],
        qc_pass_fasta = config["query_fasta"],
        designated = rules.hash_sequence_assign.output.designated,
        cached = rules.cache_sequence_assign.output.cached
    output:
        csv= os.path.join(config["tempdir"],"pangolearn_assignments.csv")
    run:

        with open(output[0],"w") as fw:
            fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
            passed = []

            version = f"PANGO-{config['pango_version']}"
            with open(input.designated,"r") as f:
                reader = csv.DictReader(f)
                note = "Assigned from designation hash."
                for row in reader:
                    
                    fw.write(f"{row['taxon']},{row['lineage']},,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},passed_qc,{note}\n")
                    passed.append(row['taxon'])

            with open(input.cached,"r") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    fw.write(f"{row['taxon']},{row['lineage']},{row['conflict']},{row['ambiguity_score']},{row['scorpio_call']},{row['scorpio_support']},{row['scorpio_conflict']},{row['version']},{row['pangolin_version']},{row['pangoLEARN_version']},{row['pango_version']},{row['status']},{row['note']}\n")
                    passed.append(row['taxon'])

            version = f"PLEARN-{config['pango_version']}"
            with open(input.qcpass, "r") as f:
                reader = csv.DictReader(f)

                for row in reader:
                    note = ''

                    support =  1 - round(float(row["score"]), 2)
                    
                    non_zero_ids = row["non_zero_ids"].split(";")
                    if len(non_zero_ids) > 1:
                        note = f"Alt assignments: {row['non_zero_ids']},{row['non_zero_scores']}"
                    
                    fw.write(f"{row['taxon']},{row['prediction']},{support},{row['imputation_score']},,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},passed_qc,{note}\n")
                    passed.append(row['taxon'])
            
            version = f"PANGO-{config['pango_version']}"
            for record in SeqIO.parse(input.qcfail,"fasta"):
                desc_list = record.description.split(" ")
                note = ""
                for i in desc_list:
                    if i.startswith("fail="):
                        note = i.lstrip("fail=")

                fw.write(f"{record.id},None,,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},fail,{note}\n")
            
            for record in SeqIO.parse(input.qc_pass_fasta,"fasta"):
                if record.id not in passed:
                    fw.write(f"{record.id},{UNASSIGNED_LINEAGE_REPORTED},,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},fail,failed_to_map\n")

rule generate_report:
    input:
        csv = os.path.join(config["tempdir"],"pangolearn_assignments.csv"),
        scorpio_voc_report = rules.scorpio.output.report,
        constellations_list = rules.get_constellations.output.list,
        alias_file = config["alias_file"]
    output:
        csv = config["outfile"]
    run:
        voc_list = []
        with open(input.constellations_list,"r") as f:
            for line in f:
                voc_list.append(line.rstrip())

        voc_dict = {}
        with open(input.scorpio_voc_report,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                if row["constellations"] != "":
                    voc_dict[row["query"]] = row

        alias_dict = {}
        with open(input.alias_file, "r") as read_file:
            alias_dict = json.load(read_file)
        if "A" in alias_dict:
            del alias_dict["A"]
        if "B" in alias_dict:
            del alias_dict["B"]

        with open(output.csv, "w") as fw:

            with open(input.csv, "r") as f:
                reader = csv.DictReader(f)
                header_names = reader.fieldnames
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()

                for row in reader:
                    new_row = row
                    if row["taxon"] in voc_dict:
                        scorpio_call_info = voc_dict[row["taxon"]]
                        new_row["scorpio_call"] = scorpio_call_info["constellations"]
                        new_row["scorpio_support"] = scorpio_call_info["support"]
                        new_row["scorpio_conflict"] = scorpio_call_info["conflict"]
                        new_row["note"] = f'scorpio call: Alt alleles {scorpio_call_info["alt_count"]}; Ref alleles {scorpio_call_info["ref_count"]}; Amb alleles {scorpio_call_info["ambig_count"]}; Oth alleles {scorpio_call_info["other_count"]}'

                        scorpio_lineage = scorpio_call_info["mrca_lineage"]
                        expanded_scorpio_lineage = expand_alias(scorpio_lineage, alias_dict)
                        expanded_pango_lineage = expand_alias(row['lineage'], alias_dict)
                        if '/' not in scorpio_lineage:
                            if expanded_scorpio_lineage and expanded_pango_lineage and not expanded_pango_lineage.startswith(expanded_scorpio_lineage):
                                new_row["note"] += f'; scorpio replaced lineage assignment {row["lineage"]}'
                                new_row['lineage'] = scorpio_lineage
                            elif "incompatible_lineages" in scorpio_call_info and row['lineage'] in scorpio_call_info["incompatible_lineages"].split("|"):
                                new_row["note"] += f'; scorpio replaced lineage assignment {row["lineage"]}'
                                new_row['lineage'] = scorpio_lineage
                    else:
                        expanded_pango_lineage = expand_alias(row['lineage'], alias_dict)
                        while expanded_pango_lineage and len(expanded_pango_lineage) > 3:
                            for voc in voc_list:
                                if expanded_pango_lineage.startswith(voc + ".") or expanded_pango_lineage == voc:
                                    # have no scorpio call but a pangolearn voc/vui call
                                    new_row['note'] += f'pangoLEARN lineage assignment {row["lineage"]} was not supported by scorpio'
                                    new_row['lineage'] = UNASSIGNED_LINEAGE_REPORTED
                                    new_row['conflict'] = ""
                                    new_row['ambiguity_score'] = ""
                                    break
                            if new_row['lineage'] == UNASSIGNED_LINEAGE_REPORTED:
                                break
                            expanded_pango_lineage = ".".join(expanded_pango_lineage.split(".")[:-1])
                    writer.writerow(new_row)

        print(green(f"Output file written to: ") + f"{output.csv}")
        if config["alignment_out"]:
            print(green(f"Output alignment written to: ") + config["outdir"] +"/sequences.aln.fasta")
