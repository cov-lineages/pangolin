#!/usr/bin/env python

import re

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


def generate_final_report():
    alias_dict = get_alias_dict(input.alias_file)

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



def scorpio_parsing():
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
        elif "incompatible_lineages" in scorpio_call_info and lineage in scorpio_call_info["incompatible_lineages"].split("|"):
            note += f'; scorpio replaced lineage assignment {lineage}'
            lineage = scorpio_lineage

def usher_parsing(usher_result,output_report):
    """
    Parsing the output of usher inference into a lineage report with columns
    name, lineage, conflict, histogram_note
    """

    with open(output_report, "w") as fw:
        fw.write("taxon,lineage,conflict,histogram_note\n")

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