#!/usr/bin/env python

import os
import re
import csv
import json

from pangolin.utils.config import *

def usher_parsing(usher_result,output_report):
    """
    Parsing the output of usher inference into a lineage report with columns
    hash, lineage, conflict, usher_note
    """

    with open(output_report, "w") as fw:
        fw.write("hash,lineage,conflict,usher_note\n")

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
    hash, lineage, conflict, ambiguity_score, pangolearn_note
    """
    with open(output_report,"w") as fw:
        fw.write("hash,lineage,conflict,ambiguity_score,pangolearn_note\n")

        with open(pangolearn_result, "r") as f:
            reader = csv.DictReader(f)

            for row in reader:
                note = ''
                support =  round((1 - float(row["score"])), 2)
                ambiguity_score = round(float(row['imputation_score']), 2)

                non_zero_ids = row["non_zero_ids"].split(";")
                if len(non_zero_ids) > 1:
                    note = f"Alt assignments: {row['non_zero_ids']},{row['non_zero_scores']}"
                
                fw.write(f"{row['taxon']},{row['prediction']},{support},{ambiguity_score},{note}\n")


def expand_alias(pango_lineage, alias_dict):
    if not pango_lineage or pango_lineage in ["None", None, "", UNASSIGNED_LINEAGE_REPORTED] or "/" in pango_lineage:
        return None

    lineage_parts = pango_lineage.split(".")
    while lineage_parts[0] in alias_dict.keys() and not lineage_parts[0].startswith('X'):
        if len(lineage_parts) > 1:
            pango_lineage = alias_dict[lineage_parts[0]] + "." + ".".join(lineage_parts[1:])
        else:
            pango_lineage = alias_dict[lineage_parts[0]]
        lineage_parts = pango_lineage.split(".")
    if lineage_parts[0] not in ["A","B"] and not lineage_parts[0].startswith('X'):
        return None
    return pango_lineage

def get_recombinant_parents(pango_lineage, alias_dict):
    # return the direct parents and the expanded lineage
    list_expanded = []
    if pango_lineage in ["None", None, "", UNASSIGNED_LINEAGE_REPORTED]:
        return list_expanded

    lineage_parts = pango_lineage.split(".")
    if not lineage_parts[0].startswith('X'):
        return list_expanded

    for parent in alias_dict[lineage_parts[0]]:
        list_expanded.append(expand_alias(parent, alias_dict))

    return list_expanded

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
            inference_dict[row["hash"]] = row
    return inference_dict

def get_cached_dict(cached_csv):
    cached_dict = {}
    if os.path.exists(cached_csv):
        with open(cached_csv, "r") as f:
            reader = csv.DictReader(f)
            for row in reader:
                cached_dict[row["hash"]] = row
    return cached_dict

def get_voc_list(voc_file, alias_file):
    voc_list = []
    alias_dict = get_alias_dict(alias_file)

    with open(voc_file,"r") as f:
        for line in f:
            expanded_voc = expand_alias(line.rstrip(), alias_dict)
            if expanded_voc:
                voc_list.append(expanded_voc)
    return voc_list

def add_relevant_fields_to_new_row(data_dict,new_row):
    for field in data_dict:
        if field in FINAL_HEADER:
            new_row[field] = data_dict[field]
        elif field in HEADER_FIELD_MAP:
            final_field = HEADER_FIELD_MAP[field]
            new_row[final_field] = data_dict[field]

def append_note(new_row, new_note):
    if new_row["note"]:
        new_row["note"] += "; " + new_note
    else:
        new_row["note"] = new_note

def generate_final_report(preprocessing_csv, inference_csv, cached_csv, alias_file, voc_list, pango_version, analysis_mode, skip_cache, output_report, config):
    """
    preprocessing_csv header is: 
    ["name","hash","lineage","scorpio_constellations",
    "scorpio_mrca_lineage","scorpio_incompatible_lineages",
    "scorpio_support","scorpio_conflict","scorpio_notes",
    "designated","qc_status","qc_notes"]

    inference_csv header is:
    usher: hash,lineage,conflict,usher_note
    pangolearn: hash,lineage,conflict,ambiguity_score,pangolearn_note

    cached_csv header is:
    hash,lineage,conflict,version,note
    """
    # the lineage aliases
    alias_dict = get_alias_dict(alias_file)

    # the output(s) of pangolearn/usher inference pipelines
    # only pass qc records present in this file
    inference_dict = get_inference_dict(inference_csv)
    cached_dict = get_cached_dict(cached_csv)

    if analysis_mode == "pangolearn":
        version = f"PLEARN-v{pango_version}"
    elif analysis_mode == "usher":
        version = f"PUSHER-v{pango_version}"

    with open(output_report, "w") as fw:
        # the output of preprocessing csv, all records present in this file
        with open(preprocessing_csv, "r") as f:
            reader = csv.DictReader(f)

            out_header = FINAL_HEADER
            if config['expanded_lineage']:
                out_header.append('expanded_lineage')
            
            writer = csv.DictWriter(fw, fieldnames=out_header, lineterminator="\n")
            writer.writeheader()

            for row in reader:
                new_row = {}
                add_relevant_fields_to_new_row(row,new_row)
                has_assignment = False

                inference_lineage = UNASSIGNED_LINEAGE_REPORTED
                if row["hash"] in cached_dict:
                    has_assignment = True
                    cached_out = cached_dict[row["hash"]]
                    add_relevant_fields_to_new_row(cached_out, new_row)
                    inference_lineage = cached_out["lineage"]
                else:
                    new_row["version"] = version
                    if row["hash"] in inference_dict:
                        has_assignment = True
                        inference_out = inference_dict[row["hash"]]
                        add_relevant_fields_to_new_row(inference_out,new_row)
                        inference_lineage = inference_out["lineage"]

                new_row[KEY_PANGOLIN_VERSION] = config[KEY_PANGOLIN_VERSION]
                new_row[KEY_SCORPIO_VERSION] = config[KEY_SCORPIO_VERSION]
                new_row[KEY_CONSTELLATIONS_VERSION] = config[KEY_CONSTELLATIONS_VERSION]

                # if it passed qc and mapped
                if has_assignment:
                    expanded_pango_lineage = expand_alias(inference_lineage, alias_dict)

                    #1. check if hash assigned
                    if row["designated"] == "True" and not skip_cache:
                        new_row["note"] = "Assigned from designation hash."
                        new_row["version"] = f"PANGO-v{pango_version}"
                        new_row["lineage"] = row["lineage"] # revert back to designation hash lineage

                    #2. check if scorpio assigned
                    elif row["scorpio_constellations"]:
                        scorpio_lineage = row["scorpio_mrca_lineage"]
                        expanded_scorpio_lineage = expand_alias(scorpio_lineage, alias_dict)

                        recombinant_parents = get_recombinant_parents(expanded_pango_lineage, alias_dict)

                        if '/' not in scorpio_lineage:
                            if expanded_scorpio_lineage and \
                                    not expanded_pango_lineage.startswith(expanded_scorpio_lineage) and \
                                    not (expanded_pango_lineage.startswith("X") and expanded_scorpio_lineage in recombinant_parents):
                                append_note(new_row, f'scorpio replaced lineage inference {inference_lineage}')
                                new_row["lineage"] = scorpio_lineage

                            elif row["scorpio_incompatible_lineages"] and inference_lineage in row["scorpio_incompatible_lineages"].split("|"):
                                append_note(new_row, f'scorpio replaced lineage inference {inference_lineage}')
                                new_row["lineage"] = scorpio_lineage

                            elif not expanded_scorpio_lineage:
                                append_note(new_row, f'scorpio replaced lineage inference {inference_lineage}')
                                new_row['lineage'] = UNASSIGNED_LINEAGE_REPORTED

                    #3. check if lineage is a voc
                    elif row["lineage"] in voc_list:
                        while expanded_pango_lineage and len(expanded_pango_lineage) > 3:
                            for voc in voc_list:
                                if expanded_pango_lineage.startswith(voc + ".") or expanded_pango_lineage == voc:
                                    # have no scorpio call but an inference voc/vui call
                                    append_note(new_row, f'Lineage inference {inference_lineage} was not supported by scorpio')
                                    new_row['lineage'] = UNASSIGNED_LINEAGE_REPORTED
                                    new_row['conflict'] = ""
                                    new_row['ambiguity_score'] = ""
                                    break

                            if new_row['lineage'] == UNASSIGNED_LINEAGE_REPORTED:
                                break

                            expanded_pango_lineage = ".".join(expanded_pango_lineage.split(".")[:-1])

                else:
                    new_row["lineage"] = UNASSIGNED_LINEAGE_REPORTED
                if config['expanded_lineage']:
                    if new_row["lineage"] == UNASSIGNED_LINEAGE_REPORTED:
                        new_row["expanded_lineage"] = UNASSIGNED_LINEAGE_REPORTED
                    else:
                        new_row["expanded_lineage"] = expand_alias(new_row["lineage"], alias_dict)
                writer.writerow(new_row)
