#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
import pangofunks as pfunk

##### Configuration #####

if config.get("trained_model"):
    config["trained_model"] = os.path.join(workflow.current_basedir,'..', config["trained_model"])

if config.get("header_file"):
    config["header_file"] = os.path.join(workflow.current_basedir,'..', config["header_file"])

##### Target rules #####

if config.get("lineages_csv"):
    print("Going to run the global report summary")
else:
    config["lineages_csv"]=""

if config.get("usher_protobuf"):
    print("Using UShER as inference engine with tree " + config["usher_protobuf"])
else:
    config["usher_protobuf"]=""
if not config.get("threads"):
    config["threads"]=""

ruleorder: use_usher > overwrite

if config["lineages_csv"] != "":
    rule all:
        input:
            config["outfile"],
            os.path.join(config["outdir"],"global_lineage_information.csv"),
            os.path.join(config["tempdir"],"VOC_report.scorpio.csv")
else:
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
        trim_end = 29674
    output:
        fasta = os.path.join(config["aligndir"],"sequences.aln.fasta")
    log:
        os.path.join(config["tempdir"], "logs/minimap2_sam.log")
    shell:
        """
        minimap2 -a -x asm5 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} '{input.fasta}' | \
        gofasta sam toMultiAlign \
            -t {workflow.cores} \
            --reference {input.reference:q} \
            --trimstart {params.trim_start} \
            --trimend {params.trim_end} \
            --trim \
            --pad > '{output.fasta}'
        """

rule pangolearn:
    input:
        fasta = rules.align_to_reference.output.fasta,
        model = config["trained_model"],
        header = config["header_file"],
        reference = config["reference_fasta"]
    output:
        os.path.join(config["tempdir"],"lineage_report.pass_qc.csv")
    shell:
        # should output a csv file with no headers but with columns similar to:
        # f.write("taxon,prediction,score,imputation_score,non_zero_ids,non_zero_scores,designated\n")
        """
        pangolearn.py --header-file {input.header:q} --model-file {input.model:q} --reference-file {input.reference:q} --fasta {input.fasta:q} -o {output[0]:q}
        """

rule add_failed_seqs:
    input:
        qcpass= os.path.join(config["tempdir"],"lineage_report.pass_qc.csv"),
        qcfail= config["qc_fail"],
        qc_pass_fasta = config["query_fasta"]
    params:
        pangolearn_version = config["pangoLEARN_version"],
        designation_version = config["pango_version"],
        pangolin_version = config["pangolin_version"]
    output:
        csv= os.path.join(config["tempdir"],"pangolearn_assignments.csv")
    run:

        fw = open(output[0],"w")
        fw.write("taxon,lineage,voc_scorpio,conflict,ambiguity_score,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
        passed = []
        version = f"PLEARN-{params.designation_version}"
        with open(input.qcpass, "r") as f:
            reader = csv.DictReader(f)

            for row in reader:
                note = ''

                support =  1 - round(float(row["score"]), 2)
                
                non_zero_ids = row["non_zero_ids"].split(";")
                if len(non_zero_ids) > 1:
                    note = f"Alt assignments: {row['non_zero_ids']},{row['non_zero_scores']}"
                
                if row["designated"]:
                    version = f"PANGO-{params.designation_version}"
                else:
                    version = f"PLEARN-{params.designation_version}"

                fw.write(f"{row['taxon']},{row['prediction']},,{support},{row['imputation_score']},{version},{params.pangolin_version},{params.pangolearn_version},{params.designation_version},passed_qc,{note}\n")
                passed.append(row['taxon'])

        for record in SeqIO.parse(input.qcfail,"fasta"):
            desc_list = record.description.split(" ")
            note = ""
            for i in desc_list:
                if i.startswith("fail="):
                    note = i.lstrip("fail=")
            # needs to mirror the structure of the output from pangolearn
            
            #"taxon,lineage,conflict,ambiguity_score,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n"
            fw.write(f"{record.id},None,,NA,NA,{version},{params.pangolin_version},{params.pangolearn_version},{params.designation_version},fail,{note}\n")
        


        for record in SeqIO.parse(input.qc_pass_fasta,"fasta"):
            if record.id not in passed:
                fw.write(f"{record.id},None,,NA,NA,{version},{params.pangolin_version},{params.pangolearn_version},{params.designation_version},fail,failed_to_map\n")

        fw.close()

rule scorpio:
    input:
        fasta = rules.align_to_reference.output.fasta,
    output:
        report = os.path.join(config["tempdir"],"VOC_report.scorpio.csv")
    shell:
        """
        scorpio classify -i {input.fasta:q} -o {output.report:q} 
        """

rule overwrite:
    input:
        csv = os.path.join(config["tempdir"],"pangolearn_assignments.csv"),
        scorpio_voc_report = rules.scorpio.output.report
    output:
        csv = config["outfile"]
    run:
        # b117 = {}
        # with open(input.b117_variants, "r") as f:
        #     reader = csv.DictReader(f)
        #     for row in reader:
        #         if int(row["alt_count"]) > 4 and int(row["ref_count"])<6:
        #             b117[row["query"]] = {'alt':row["alt_count"], 'ref':row["ref_count"], 'oth':row["other_count"]}
        # b1351 = {}
        # with open(input.b1351_variants, "r") as f:
        #     reader = csv.DictReader(f)
        #     for row in reader:
        #         if int(row["alt_count"]) > 4 and int(row["ref_count"])<2:
        #             b1351[row["query"]] = {'alt':row["alt_count"], 'ref':row["ref_count"], 'oth':row["other_count"]}
        # p1 = {}
        # with open(input.p1_variants, "r") as f:
        #     reader = csv.DictReader(f)
        #     for row in reader:
        #         if int(row["alt_count"]) > 10:
        #             p1[row["query"]] = {'alt':row["alt_count"], 'ref':row["ref_count"], 'oth':row["other_count"]}
        # p2 = {}
        # with open(input.p2_variants, "r") as f:
        #     reader = csv.DictReader(f)
        #     for row in reader:
        #         if int(row["alt_count"]) > 4 and int(row["ref_count"])<4:
        #             p2[row["query"]] = {'alt':row["alt_count"], 'ref':row["ref_count"], 'oth':row["other_count"]}
        # p3 = {}
        # with open(input.p3_variants, "r") as f:
        #     reader = csv.DictReader(f)
        #     for row in reader:
        #         if int(row["alt_count"]) > 8 and int(row["ref_count"])<4:
        #             p3[row["query"]] = {'alt':row["alt_count"], 'ref':row["ref_count"], 'oth':row["other_count"]}
        voc_dict = {}
        with open(input.scorpio_voc_report,"r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                voc_dict[row["query"]] = row["constellations"]

        with open(output.csv, "w") as fw:
            # "taxon,lineage,probability,pangoLEARN_version,status,note" 
            with open(input.csv, "r") as f:
                reader = csv.DictReader(f)
                header_names = reader.fieldnames
                writer = csv.DictWriter(fw, fieldnames=header_names,lineterminator='\n')
                writer.writeheader()

                for row in reader:
                    new_row = row
                    if row["taxon"] in voc_dict:
                        new_row["voc_scorpio"] = voc_dict[row["taxon"]]
                    writer.writerow(new_row)

                    # if row["lineage"] =="B.1.1.7" and row["taxon"] not in b117:
                    #     new_row = row
                        
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1.1"

                    #     writer.writerow(new_row)

                    # elif row["taxon"] in b117:
                    #     new_row = row
                        
                    #     snps = b117[row["taxon"]]
                    #     note = f'{snps["alt"]}/17 B.1.1.7 SNPs ({snps["ref"]} ref and {snps["oth"]} other)'

                    #     new_row["note"] = note
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1.1.7"

                    #     writer.writerow(new_row)
                    # elif row["lineage"].startswith("B.1.351") and row["taxon"] not in b1351:
                    #     new_row = row
                        
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1"

                    #     writer.writerow(new_row)
                        
                    # elif row["taxon"] in b1351:
                    #     new_row = row
                        
                    #     snps = b1351[row["taxon"]]
                    #     note = f'{snps["alt"]}/9 B.1.351 SNPs ({snps["ref"]} ref and {snps["oth"]} other)'

                    #     new_row["note"] = note
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1.351"

                    #     writer.writerow(new_row)
                    # elif row["taxon"] in p2:
                    #     new_row = row
                        
                    #     snps = p2[row["taxon"]]
                    #     note = f'{snps["alt"]}/5 P.2 (B.1.1.28.2) SNPs ({snps["ref"]} ref and {snps["oth"]} other)'

                    #     new_row["note"] = note
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "P.2"

                    #     writer.writerow(new_row)
                    # elif row["lineage"] =="P.2" and row["taxon"] not in p2:
                    #     new_row = row
                        
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1.1.28"

                    #     writer.writerow(new_row)
                    # elif row["taxon"] in p1:
                    #     new_row = row
                        
                    #     snps = p1[row["taxon"]]
                    #     note = f'{snps["alt"]}/16 P.1 (B.1.1.28.1) SNPs ({snps["ref"]} ref and {snps["oth"]} other)'

                    #     new_row["note"] = note
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "P.1"

                    #     writer.writerow(new_row)
                    # elif row["lineage"] =="P.1" and row["taxon"] not in p1:
                    #     new_row = row
                        
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1.1.28"

                    #     writer.writerow(new_row)
                    # elif row["taxon"] in p3:
                    #     new_row = row
                    #     snps = p3[row["taxon"]]
                    #     note = f'{snps["alt"]}/12 P.3 (B.1.1.28.3) SNPs ({snps["ref"]} ref and {snps["oth"]} other)'

                    #     new_row["note"] = note
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "P.3"

                    #     writer.writerow(new_row)
                    # elif row["lineage"] =="P.3" and row["taxon"] not in p3:
                    #     new_row = row
                        
                    #     new_row["conflict"] = "0"
                    #     new_row["lineage"] = "B.1.1.28"

                    #     writer.writerow(new_row)
                    # else:
                    #     writer.writerow(row)
        print(pfunk.green(f"Output file written to: ") + f"{output.csv}")
        if config["alignment_out"]:
            print(pfunk.green(f"Output alignment written to: ") + config["outdir"] +"/sequences.aln.fasta")
rule report_results:
    input:
        csv = config["outfile"],
        lineages_csv = config["lineages_csv"]
    output:
        os.path.join(config["outdir"],"global_lineage_information.csv")
    shell:
        """
        report_results.py \
        -p {input.csv:q} \
        -b {input.lineages_csv:q} \
        -o {output:q} 
        """

rule use_usher:
    input:
        fasta = rules.align_to_reference.output.fasta,
        reference = config["reference_fasta"],
        usher_protobuf = config["usher_protobuf"]
    params:
        designation_version = config["pango_version"],
        pangolin_version = config["pangolin_version"],
        tempdir = config["tempdir"],
        threads = config["threads"],
        version = config["pangoLEARN_version"],
        vcf = os.path.join(config["tempdir"], "sequences.aln.vcf"),
    output:
        csv = config["outfile"]
    shell:
        """
        faToVcf <(cat {input.reference:q} <(echo "") {input.fasta:q}) {params.vcf}
        if [ -z "{params.threads}" ]; then
            T=""
        else
            T="-T {params.threads}"
        fi
        echo usher -i {input.usher_protobuf:q} -v {params.vcf} $T -d {params.tempdir}
        usher -i {input.usher_protobuf:q} -v {params.vcf} $T -d {params.tempdir} &> usher.log
        echo "taxon,lineage,conflict,ambiguity_score,version,pangolin_version,pangoLEARN_version,pango_version,status,note" > {output.csv}
        awk -F'\t' -v OFS=, '{{ print $1, $NF, "NA", "NA", "PUSHER-{params.designation_version}", "{params.pangolin_version}", "{params.version}", "{params.designation_version}", "passed_qc", "UShER inference"; }}' \
            clades.txt >> {output.csv}
        echo ""
        echo "Output written to {output.csv}"
        """
#"taxon,lineage,conflict,ambiguity_score,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n"