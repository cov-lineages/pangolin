#!/usr/bin/env python3
import pkg_resources
from pangolin.utils.log_colours import green,cyan,red
import sys
import os
import gzip

from pangolin.utils.config import *

def package_data_check(filename,directory,key,config):
    try:
        package_datafile = os.path.join(directory,filename)
        data = pkg_resources.resource_filename('pangolin', package_datafile)
        config[key] = data
    except:
        sys.stderr.write(cyan(f'Error: Missing package data.')+f'\n\t- {filename}\nPlease install the latest pangolin version with `pangolin --update`.\n')
        sys.exit(-1)


def check_install(config):
    resources = [
        {"key":"reference_fasta",
        "directory":"data",
        "filename":"reference.fasta"}
    ]
    for resource in resources:
        package_data_check(resource["filename"],resource["directory"],resource["key"],config)

def find_designation_cache_and_alias(datadir,designation_cache_file,alias_file):
    designation_cache = ""
    alias = ""
    for r,d,f in os.walk(datadir):
        for fn in f:
            if fn == designation_cache_file:
                designation_cache = os.path.join(r, fn)
            elif fn == alias_file:
                alias = os.path.join(r, fn)
    if designation_cache == "":
        sys.stderr.write(cyan(f'Error: Missing designation cache file. Either supply a datadir with a {designation_cache_file} file, or specify `--skip-designation-cache`\n'))
        sys.exit(-1)
    elif alias == "":
        sys.stderr.write(cyan(f'Error: Missing alias file. Please supply a datadir with a {alias_file} file or check installation of pangolin-data dependency.\n'))
        sys.exit(-1)
    return designation_cache,alias

def check_file_arg(arg_file, cwd, description):
    if arg_file:
        file_path = os.path.join(cwd, arg_file)
        if not os.path.exists(file_path):
            sys.stderr.write(f"Error: cannot find {description} file at {file_path}\n")
            sys.exit(-1)
    return file_path

def get_datafiles(datadir,file_dict,config):
    datafiles = {}
    for r,d,f in os.walk(datadir):
        for fn in f:
            if fn in file_dict:
                datafiles[file_dict[fn]] = os.path.join(r, fn)
    for fn in datafiles:
        config[fn] = datafiles[fn]
    for fn in file_dict:
        if file_dict[fn] not in config:
            sys.stderr.write(cyan(f'Error: Cannot find {fn} in datadir. Please supply a datadir with required files or specify an alternative analysis mode.\nPlease see https://cov-lineages.org/pangolin.html for full installation and updating instructions.'))
            sys.exit(-1)

    print(green("****\nData files found:"))
    for fn in datafiles:
        print(f"{fn}:\t{datafiles[fn]}")
        config[fn] = datafiles[fn]
    print(green("****"))


def install_error(package, url):
    sys.stderr.write(cyan(f'Error: please install `{package}` with \n') +
                     f"pip install git+{url}\n")
    sys.exit(-1)


def get_assignment_cache(cache_file, config):
    cache = ""
    if config[KEY_PANGOLIN_ASSIGNMENT_VERSION] is not None:
        pangolin_assignment_dir = config[KEY_PANGOLIN_ASSIGNMENT_PATH]
        for r, d, f in os.walk(pangolin_assignment_dir):
            for fn in f:
                if fn == cache_file and cache == "":
                    cache = os.path.join(r, fn)
        if not os.path.exists(cache):
            sys.stderr.write(cyan(f'Error: cannot find assignment cache file {cache_file} in pangolin_assignment\n'))
            sys.exit(-1)
    else:
        sys.stderr.write(cyan('\nError: "pangolin --add-assignment-cache" is required before '
                              '"pangolin --use-assignment-cache", in order to install optional '
                              'pangolin-assignment repository (that will make future data updates slower).\n'))
        sys.exit(-1)

    # Check versions of pangolin-data and pangolin-assignment to make sure they are consistent.
    if config[KEY_PANGOLIN_ASSIGNMENT_VERSION].lstrip('v') != config[KEY_PANGOLIN_DATA_VERSION].lstrip('v'):
        print(cyan(f'Error: pangolin_assignment cache version {config[KEY_PANGOLIN_ASSIGNMENT_VERSION]} '
                   f'does not match pangolin_data version {config[KEY_PANGOLIN_DATA_VERSION]}. '
                   'Run "pangolin --update-data" to fetch latest versions of both.'))
        sys.exit(-1)

    try:
        with gzip.open(cache, 'rt') as f:
            line = f.readline()
    except:
        with open(cache, 'r') as f:
            line = f.readline()
            if "git-lfs.github.com" in line:
                sys.stderr.write(cyan(
                    'Error: Git LFS file not pulled successfully. Please install git-lfs \nusing conda or an alternative (not pip) then re-install pangolin-assignment \nwith pip install git+https://github.com/cov-lineages/pangolin-assignment.git\n')
                )
                sys.exit(-1)
    return cache

def get_constellation_files(path):
    constellation_files = []
    if path is not None:
        for r, _, f in os.walk(path):
            for fn in f:
                if (r.endswith('/constellations') or r.endswith('/constellations/definitions')) and fn.endswith('.json'):
                    constellation_files.append(os.path.join(r, fn))
    return constellation_files

# config={}
# check_install()
