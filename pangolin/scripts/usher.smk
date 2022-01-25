#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
import gzip
from pangolin.utils.log_colours import green,cyan
import pangolin.utils.report_collation as report_collation
from pangolin.utils.config import *

##### Report options #####


##### Target rules #####

rule all:
    input:
        csv = os.path.join(config[KEY_TEMPDIR],"inference_report.csv")

rule usher_inference:
    input:
        fasta = os.path.join(config[KEY_TEMPDIR],"pass_qc.fasta"),
        reference = config[KEY_REFERENCE_FASTA],
        usher_protobuf = config[KEY_USHER_PB]
    params:
        vcf = os.path.join(config[KEY_TEMPDIR], "sequences.aln.vcf")
    threads: workflow.cores
    output:
        txt = os.path.join(config[KEY_TEMPDIR], "clades.txt")
    log:
        os.path.join(config[KEY_TEMPDIR], "logs/usher.log")
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
        txt = rules.usher_inference.output.txt,
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"inference_report.csv")
    run:
        report_collation.usher_parsing(input.txt, output.csv)
