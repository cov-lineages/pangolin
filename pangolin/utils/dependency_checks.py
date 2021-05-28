#!/usr/bin/env python3
import subprocess
import os
import sys
from pangolin.utils import log_colours as colour
import importlib

import pangolin.utils.custom_logger as custom_logger
import pangolin.utils.log_handler_handle as lh

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

def check_dependencies():

    missing = []

    dependency_list = ["gofasta","minimap2","snakemake","usher"]
    module_list = ["Bio","sklearn","pandas","joblib","pysam","pangoLEARN","constellations"]

    for dependency in dependency_list:
        check_this_dependency(dependency, missing)

    for module in module_list:
        check_module(module, missing)

    if missing:
        if len(missing)==1:
            sys.stderr.write(colour.cyan(f'Error: Missing dependency `{missing[0]}`.')+'\nPlease update your pangolin environment.\n')
            sys.exit(-1)
        else:
            dependencies = ""
            for i in missing:
                dependencies+=f"\t- {i}\n"

            sys.stderr.write(colour.cyan(f'Error: Missing dependencies.')+f'\n{dependencies}Please update your pangolin environment.\n')
            sys.exit(-1)
    else:
        print(colour.green("All dependencies satisfied."))

def set_up_verbosity(config):
    if config["verbose"]:
        config["quiet"] = False
        config["log_api"] = ""
        config["log_string"] = ""
    else:
        config["quiet"] = True
        logger = custom_logger.Logger()
        config["log_api"] = logger.log_handler

        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "
# check_dependencies()
