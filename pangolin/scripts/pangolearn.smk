#!/usr/bin/env python

import csv
from Bio import SeqIO
import os
import gzip
from pangolin.utils.log_colours import green,cyan,red

from pangolin.utils.report_collation import pangolearn_parsing
import pangolin.pangolearn.pangolearn as pangolearn
from pangolin.utils.config import *
##### Configuration #####


rule all:
    input:
        os.path.join(config[KEY_TEMPDIR],"inference_report.csv")

rule pangolearn:
    input:
        fasta = os.path.join(config[KEY_TEMPDIR],"pass_qc.fasta"),
        model = config[KEY_PLEARN_MODEL],
        header = config[KEY_PLEARN_HEADER],
        reference = config[KEY_REFERENCE_FASTA]
    output:
        csv = os.path.join(config[KEY_TEMPDIR],"pangoLEARN_report.csv")
    run:
        pangolearn.assign_lineage(input.header,input.model,input.reference,input.fasta,output.csv)

rule pangolearn_output:
    input:
        csv= rules.pangolearn.output.csv
    output:
        csv= os.path.join(config[KEY_TEMPDIR],"inference_report.csv")
    run:
        pangolearn_parsing(input.csv, output.csv)

