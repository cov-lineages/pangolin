
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
import subprocess

import importlib

END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'


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

def which(dependency):
    try:
        subprocess.check_output(["which", dependency])
        return True
    except subprocess.CalledProcessError:
        return False

def check_module(module, missing):
    try:
        importlib.import_module(module)
    except ImportError:
        missing.append(module)

def check_this_dependency(dependency,missing):
    check = which(dependency)

    if not check:
        missing.append(dependency)

def check_installs():

    missing = []

    dependency_list = ["gofasta","minimap2","snakemake"]
    module_list = ["Bio","sklearn","pandas","joblib","pysam","pangoLEARN"]

    for dependency in dependency_list:
        check_this_dependency(dependency, missing)

    for module in module_list:
        check_module(module, missing)

    if missing:
        if len(missing)==1:
            sys.stderr.write(cyan(f'Error: Missing dependency `{missing[0]}`.')+'\nPlease update your pangolin environment.\n')
            sys.exit(-1)
        else:
            dependencies = ""
            for i in missing:
                dependencies+=f"\t- {i}\n"

            sys.stderr.write(cyan(f'Error: Missing dependencies.')+f'\n{dependencies}Please update your pangolin environment.\n')
            sys.exit(-1)
    else:
        print(green("All dependencies satisfied."))


# check_installs()
