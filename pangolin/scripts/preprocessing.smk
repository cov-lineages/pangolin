#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
from pangolin.utils.log_colours import green,cyan
from pangolin.utils.hash_functions import get_hash_string

from pangolin.utils.preprocessing import *

from pangolin.utils.config import *

rule all:
    input:
        config[KEY_ALIGNMENT_FILE],
         os.path.join(config[KEY_TEMPDIR], "get_constellations.txt"),
        os.path.join(config[KEY_TEMPDIR],"preprocessing.csv")

rule align_to_reference:
    input:
        fasta = config[KEY_QUERY_FASTA],
        reference = config[KEY_REFERENCE_FASTA]
    params:
        trim_start = config[KEY_TRIM_START],
        trim_end = config[KEY_TRIM_END],
        sam = os.path.join(config[KEY_TEMPDIR],"mapped.sam")
    output:
        fasta = config[KEY_ALIGNMENT_FILE]
    log:
        os.path.join(config[KEY_TEMPDIR], "logs/minimap2_sam.log")
    shell:
    # the first line of this streams through the fasta and replaces '-' in sequences with empty strings
    # this could be replaced by a python script later
        """
        awk '{{ if ($0 !~ /^>/) {{ gsub("-", "",$0); }} print $0; }}' {input.fasta} | \
        minimap2 -a -x asm20 --sam-hit-only --secondary=no -t  {workflow.cores} {input.reference:q} - -o {params.sam:q} &> {log:q} 
        gofasta sam toMultiAlign \
            -s {params.sam:q} \
            -t {workflow.cores} \
            --reference {input.reference:q} \
            --trimstart {params.trim_start} \
            --trimend {params.trim_end} \
            --trim \
            --pad > '{output.fasta}'
        """

rule create_seq_hash:
    input:
        rules.align_to_reference.output.fasta
    output:
        csv =  os.path.join(config[KEY_TEMPDIR],"hash_map.csv"),
        fasta = os.path.join(config[KEY_TEMPDIR],"hashed.aln.fasta")
    run:
        record_count,collapsed_count = create_seq_hash(input[0],output.csv,output.fasta)
        print(green("****\nQuery sequences collapsed from ") + f"{record_count}" +green(" to ") + f"{collapsed_count}" + green(" unique sequences."))

rule cache_sequence_assign:
    input:
        rules.create_seq_hash.output.csv
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"designation_status.csv")
    run:
        designation_count = designation_assign(config[KEY_DESIGNATION_CACHE],input[0],output.csv)
        print(green("****\n") + f"{designation_count}" +green(" sequences assigned via designations."))

rule scorpio:
    input:
        fasta = rules.create_seq_hash.output.fasta,
    params:
        constellation_files = " ".join(config[KEY_CONSTELLATION_FILES])
    output:
        report = os.path.join(config[KEY_TEMPDIR],"VOC_report.scorpio.csv")
    threads:
        workflow.cores
    log:
        os.path.join(config[KEY_TEMPDIR], "logs/scorpio.log")
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
        constellation_files = " ".join(config[KEY_CONSTELLATION_FILES])
    output:
        list = os.path.join(config[KEY_TEMPDIR], "get_constellations.txt")
    shell:
        """
        scorpio list \
        --constellations {params.constellation_files} \
        --pangolin > {output.list:q}
        """

rule sequence_qc:
    input:
        rules.create_seq_hash.output.fasta
    output:
        pass_qc = os.path.join(config[KEY_TEMPDIR],"pass_qc.fasta"),
        csv = os.path.join(config[KEY_TEMPDIR],"seq_status.csv")
    run:
        print(green("****\nRunning sequence QC"))
        total_pass = seq_qc(input[0],output.pass_qc,output.csv,config[KEY_MAXAMBIG])
        print(green("Total passing QC: ") + f"{total_pass}")
        print(green("****"))
        # if total_pass == 0:
        #     with open(outfile, "w") as fw:
        #         fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
        #         for record in SeqIO.parse(os.path.join(tempdir,'query.failed_qc.fasta'), "fasta"):
        #             desc = record.description.split(" ")
        #             reason = ""
        #             for item in desc:
        #                 if item.startswith("fail="):
        #                     reason = item.split("=")[1]
        #             fw.write(f"{record.id},None,,,,,,PANGO-{PANGO_VERSION},{__version__},{pangoLEARN.__version__},{PANGO_VERSION},fail,{reason}\n")
        #     print(cyan(f'Note: no query sequences have passed the qc\n'))
        #     sys.exit(0)

rule merged_info:
    input:
        fasta = config[KEY_QUERY_FASTA],
        qc_status = rules.sequence_qc.output.csv,
        voc_report = rules.scorpio.output.report,
        designated = rules.cache_sequence_assign.output.csv,
        hash_map = rules.create_seq_hash.output.csv
    output:
        merged = os.path.join(config[KEY_TEMPDIR],"preprocessing.csv")
    run:
        merge_files(input.fasta,input.qc_status,input.voc_report, input.designated, input.hash_map, output.merged)

