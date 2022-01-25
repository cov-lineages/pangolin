#!/usr/bin/env python

import re
import csv

def usher_parsing(usher_result,output_report):
    """
    Parsing the output of usher inference into a lineage report with columns
    name, lineage, conflict, usher_note
    """

    with open(output_report, "w") as fw:
        fw.write("taxon,lineage,conflict,usher_note\n")

        with open(usher_result, "r") as f:
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
                
                fw.write(f"{name},{lineage},{conflict},{histogram_note}\n")

def pangolearn_parsing(pangolearn_result,output_report):
    """
    Parsing the output of pangolearn inference into a lineage report with columns
    name, lineage, conflict, ambiguity_score, pangolearn_note
    """
    with open(output_report,"w") as fw:
        fw.write("taxon,lineage,conflict,ambiguity_score,pangolearn_note\n")

        with open(pangolearn_result, "r") as f:
            reader = csv.DictReader(f)

            for row in reader:
                note = ''
                support =  1 - round(float(row["score"]), 2)
                ambiguity_score = round(float(row['imputation_score']), 2)

                non_zero_ids = row["non_zero_ids"].split(";")
                if len(non_zero_ids) > 1:
                    note = f"Alt assignments: {row['non_zero_ids']},{row['non_zero_scores']}"
                
                fw.write(f"{row['taxon']},{row['prediction']},{support},{ambiguity_score},{note}\n")


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

def get_alias_dict(alias_file):
    alias_dict = {}
    with open(alias_file, "r") as read_file:
        alias_dict = json.load(read_file)
    if "A" in alias_dict:
        del alias_dict["A"]
    if "B" in alias_dict:
        del alias_dict["B"]
    return alias_dict

def get_inference_dict(inference_csv):
    inference_dict = {}
    with open(inference_csv,"r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            inference_dict[row["taxon"]] = row
    return inference_dict

def parse_scorpio_constellation(row, inference_out, voc_list):

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

def generate_final_report(preprocessing_csv, inference_csv, alias_file, output_report):
    alias_dict = get_alias_dict(alias_file)
    inference_dict = get_inference_dict(inference_csv)
    
    with open(preprocessing_csv, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:

            lineage_note = ""
            version = ""

            if row["hash"] in inference_dict:
                inference_out = inference_dict[row["hash"]]

                #1. check if hash assigned
                if row["designated"] == "True":
                    lineage_note = "Assigned from designation hash."
                    version = f"PANGO-{config['pango_version']}"
                
                #2. check if scorpio assigned
                elif row["scorpio_call"]:
                    scorpio_lineage = row["scorpio_mrca_lineage"]
                    expanded_scorpio_lineage = expand_alias(scorpio_lineage, alias_dict)
                    expanded_pango_lineage = expand_alias(inference_out["lineage"], alias_dict)





            else:
                inference_out = None
            




    ## Catching scorpio and usher output 
    with open(output_report, "w") as fw:
        fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
        
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

        version = f"PUSHER-{config['pango_version']}"

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
    if config[KEY_ALIGNMENT_OUT]:
        print(green(f"Output alignment written to: ") + config[KEY_OUTDIR] +"/sequences.aln.fasta")


