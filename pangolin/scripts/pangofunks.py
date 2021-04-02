
#!/usr/bin/env python3

import os
import argparse
import csv 
import sys
from Bio import SeqIO
from datetime import datetime 
from datetime import date
import tempfile
import pkg_resources
import yaml
import subprocess

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'


def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'dim' in text_colour:
        coloured_text = DIM
    elif 'cyan' in text_colour:
        coloured_text = 'cyan'
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text

def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING

def bold_underline(text):
    return BOLD + UNDERLINE + text + END_FORMATTING


def check_installs():
    go_fasta_check = os.system("gofasta sam -h")
    
    if not go_fasta_check == 0:
        sys.stderr.write(cyan('Error: Missing dependency `gofasta`.')+'\nPlease update your pangolin environment or install gofasta with `conda install gofasta -c bioconda`\n')
        sys.exit(-1)

    minimap2_check = os.system("minimap2 --version")
    
    if not minimap2_check == 0:
        sys.stderr.write(cyan('Error: Missing dependency `minimap2`.')+'\nPlease update your pangolin environment\n')
        sys.exit(-1)

    snakemake_check = os.system("snakemake --version")

    if not snakemake_check == 0:
        sys.stderr.write(cyan('Error: Missing dependency `snakemake-minimal`.')+'\nPlease update your pangolin environment\n')
        sys.exit(-1)