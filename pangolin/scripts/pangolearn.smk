#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
from pangolin.utils.log_colours import green,cyan,red
from pangolin.utils.hash_functions import get_hash_string
import pangolin.pangolearn.pangolearn as pangolearn

##### Configuration #####

if config.get("trained_model"):
    config["trained_model"] = os.path.join(workflow.current_basedir,'..', config["trained_model"])

if config.get("header_file"):
    config["header_file"] = os.path.join(workflow.current_basedir,'..', config["header_file"])

##### Utility functions #####

def expand_alias(pango_lineage, alias_dict):
    if not pango_lineage or pango_lineage in ["None", None, ""] or "/" in pango_lineage:
        return None

    lineage_parts = pango_lineage.split(".")
    if lineage_parts[0].startswith('X'):
        return pango_lineage
    while lineage_parts[0] in alias_dict.keys():
        if len(lineage_parts) > 1:
            pango_lineage = alias_dict[lineage_parts[0]] + "." + ".".join(lineage_parts[1:])
        else:
            pango_lineage = alias_dict[lineage_parts[0]]
        lineage_parts = pango_lineage.split(".")
    if lineage_parts[0] not in ["A","B"]:
        return None
    return pango_lineage


##### Target rules #####

if not config.get("usher_protobuf"):
    config["usher_protobuf"]=""

ruleorder: usher_to_report > generate_report

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
        minimap2 -a -x asm5 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} '{input.fasta}' -o {params.sam:q} &> {log:q} 
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
                        if hash_string in set_hash:
                            fw.write(f"{record.id},{set_hash[hash_string]}\n")
                        else:
                            fseq.write(f">{record.description}\n{record.seq}\n")

rule pangolearn:
    input:
        fasta = rules.hash_sequence_assign.output.for_inference,
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
        designated = rules.hash_sequence_assign.output.designated
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
                    fw.write(f"{record.id},None,,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},fail,failed_to_map\n")

rule scorpio:
    input:
        fasta = rules.align_to_reference.output.fasta,
    params:
        constellations_files = " ".join(config["constellations_files"])
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
        --constellations {params.constellations_files} \
        --pangolin \
        --long &> {log:q}
        """

rule generate_report:
    input:
        csv = os.path.join(config["tempdir"],"pangolearn_assignments.csv"),
        scorpio_voc_report = rules.scorpio.output.report,
        alias_file = config["alias_file"]
    output:
        csv = config["outfile"]
    run:
        
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
                        new_row["note"] = f'scorpio call: Alt alleles {scorpio_call_info["alt_count"]}; Ref alleles {scorpio_call_info["ref_count"]}; Amb alleles {scorpio_call_info["ambig_count"]}'

                        scorpio_lineage = scorpio_call_info["mrca_lineage"]
                        expanded_scorpio_lineage = expand_alias(scorpio_lineage, alias_dict)
                        expanded_pango_lineage = expand_alias(row['lineage'], alias_dict)
                        if '/' not in scorpio_lineage:
                            if expanded_scorpio_lineage and expanded_pango_lineage and not expanded_pango_lineage.startswith(expanded_scorpio_lineage):
                                new_row["note"] += f'; scorpio replaced lineage assignment {row["lineage"]}'
                                new_row['lineage'] = scorpio_lineage

                    writer.writerow(new_row)

        print(green(f"Output file written to: ") + f"{output.csv}")
        if config["alignment_out"]:
            print(green(f"Output alignment written to: ") + config["outdir"] +"/sequences.aln.fasta")

rule use_usher:
    input:
        fasta = rules.hash_sequence_assign.output.for_inference,
        reference = config["reference_fasta"],
        usher_protobuf = config["usher_protobuf"]
    params:
        vcf = os.path.join(config["tempdir"], "sequences.aln.vcf")
    threads: workflow.cores
    output:
        txt = os.path.join(config["tempdir"], "clades.txt")
    log:
        os.path.join(config["tempdir"], "logs/usher.log")
    shell:
        """
        echo "Using UShER as inference engine."
        if [ -s {input.fasta:q} ]; then
            faToVcf <(cat {input.reference:q} <(echo "") {input.fasta:q}) {params.vcf:q}
            usher -n -D -i {input.usher_protobuf:q} -v {params.vcf:q} -T {workflow.cores} -d '{config[tempdir]}' &> {log}
        else
            rm -f {output.txt:q}
            touch {output.txt:q}
        fi
        """

rule usher_to_report:
    input:
        txt = rules.use_usher.output.txt,
        scorpio_voc_report = rules.scorpio.output.report,
        designated = rules.hash_sequence_assign.output.designated,
        qcfail= config["qc_fail"],
        qc_pass_fasta = config["query_fasta"],
        alias_file = config["alias_file"]
    output:
        csv = config["outfile"]
    run:
        voc_dict = {}
        passed = []
        

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

        ## Catching scorpio and usher output 
        with open(output.csv, "w") as fw:
            fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
            
            version = f"PANGO-{config['pango_version']}"
            with open(input.designated,"r") as f:
                reader = csv.DictReader(f)
                note = "Assigned from designation hash."
                for row in reader:
                    
                    fw.write(f"{row['taxon']},{row['lineage']},,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},passed_qc,{note}\n")
                    passed.append(row['taxon'])

            version = f"PUSHER-{config['pango_version']}"
            with open(input.txt, "r") as f:
                for l in f:
                    name,lineage_histogram = l.rstrip("\n").split("\t")
                    if "*|" in lineage_histogram:
                        # example: A.28*|A.28(1/10),B.1(6/10),B.1.511(1/10),B.1.518(2/10)
                        lineage,histogram = lineage_histogram.split("*|")
                        histo_list = [ i for i in histogram.split(",") if i ]
                        conflict = 0.0
                        if len(histo_list) > 1:
                            max_count = 0
                            max_lineage = ""
                            selected_count = 0
                            total = 0
                            for lin_counts in histo_list:
                                m = re.match('([A-Z0-9.]+)\(([0-9]+)/([0-9]+)\)', lin_counts)
                                if m:
                                    lin, place_count, total = [m.group(1), int(m.group(2)), int(m.group(3))]
                                    if place_count > max_count:
                                        max_count = place_count
                                        max_lineage = lin
                                    if lin == lineage:
                                        selected_count = place_count
                            if selected_count < max_count:
                                # The selected placement was not in the lineage with the plurality
                                # of placements; go with the plurality.
                                lineage = max_lineage
                                conflict = (total - max_count) / total
                            elif total > 0:
                                conflict = (total - selected_count) / total
                        histogram_note = "Usher placements: " + " ".join(histo_list)
                    else:
                        lineage = lineage_histogram
                        conflict = ""
                        histogram_note = ""
                    scorpio_call_info,scorpio_call,scorpio_support,scorpio_conflict,note='','','','',''
                    if name in voc_dict:
                        scorpio_call_info = voc_dict[name]
                        scorpio_call = scorpio_call_info["constellations"]
                        scorpio_support = scorpio_call_info["support"]
                        scorpio_conflict = scorpio_call_info["conflict"]
                        note = f'scorpio call: Alt alleles {scorpio_call_info["alt_count"]}; Ref alleles {scorpio_call_info["ref_count"]}; Amb alleles {scorpio_call_info["ambig_count"]}'

                        scorpio_lineage = scorpio_call_info["mrca_lineage"]
                        expanded_scorpio_lineage = expand_alias(scorpio_lineage, alias_dict)
                        expanded_pango_lineage = expand_alias(lineage, alias_dict)
                        if expanded_scorpio_lineage and expanded_pango_lineage and not expanded_pango_lineage.startswith(expanded_scorpio_lineage):
                            note += f'; scorpio replaced lineage assignment {lineage}'
                            lineage = scorpio_lineage

                        if histogram_note:
                            note += f'; {histogram_note}'
                    else:
                        note = histogram_note
                    fw.write(f"{name},{lineage},{conflict},,{scorpio_call},{scorpio_support},{scorpio_conflict},{version},{config['pangolin_version']},,{config['pango_version']},passed_qc,{note}\n")
                    passed.append(name)

            version = f"PANGO-{config['pango_version']}"
            ## Catching sequences that failed qc in the report
            for record in SeqIO.parse(input.qcfail,"fasta"):
                desc_list = record.description.split(" ")
                note = ""
                for i in desc_list:
                    if i.startswith("fail="):
                        note = i.lstrip("fail=")

                fw.write(f"{record.id},None,,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},fail,{note}\n")
            
            for record in SeqIO.parse(input.qc_pass_fasta,"fasta"):
                if record.id not in passed:
                    fw.write(f"{record.id},None,,,,,,{version},{config['pangolin_version']},{config['pangoLEARN_version']},{config['pango_version']},fail,failed_to_map\n")

        print(green(f"Output file written to: ") + f"{output.csv}")
        if config["alignment_out"]:
            print(green(f"Output alignment written to: ") + config["outdir"] +"/sequences.aln.fasta")
