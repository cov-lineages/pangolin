#!/usr/bin/env python3
import pkg_resources
from pangolin.utils.log_colours import green,cyan,red
import sys
import os
import gzip

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

def find_designation_cache(datadir,designation_cache_file):
    designation_cache = ""
    for r,d,f in os.walk(datadir):
        for fn in f:
            if fn == designation_cache_file:
                designation_cache = os.path.join(r, fn)
    if designation_cache == "":
        sys.stderr.write(cyan(f'Error: Missing designation cache file. Either supply a datadir with a {designation_cache_file} file, or specify `--skip-designation-cache`\n'))
        sys.exit(-1)
    
    return designation_cache

def get_usher_protobuf_arg(usher_arg,cwd):
    if usher_arg:
        usher_protobuf = os.path.join(cwd, usher_arg)
        if not os.path.exists(usher_protobuf):
            sys.stderr.write('Error: cannot find --usher-tree file at {}\n'.format(usher_protobuf))
            sys.exit(-1)

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


def get_cache():
    cache = ""
    try:
        import pangolin_assignment
        pangolin_assignment_dir = pangolin_assignment.__path__[0]
        for r, d, f in os.walk(pangolin_assignment_dir):
            for fn in f:
                if fn == "pango_assignment.cache.csv.gz" and cache == "":
                    cache = os.path.join(r, fn)
        if not os.path.exists(cache):
            sys.stderr.write('Error: cannot find cache\n')
            sys.exit(-1)
    except:
        sys.stderr.write(cyan('Error: please install `pangolin_assignment` with \n') +
                            "pip install git+https://github.com/cov-lineages/pangolin-assignment.git")
        sys.exit(-1)

    try:
        with gzip.open(cache, 'rt') as f:
            line = f.readline()
    except:
        with open(cache, 'r') as f:
            line = f.readline()
            if "git-lfs.github.com" in line:
                sys.stderr.write(
                    'Error: Git LFS file not pulled successfully. Please install git-lfs \nusing conda or an alternative (not pip) then re-install pangolin-assignment \nwith pip install git+https://github.com/cov-lineages/pangolin-assignment.git\n')
                sys.exit(-1)
    return cache

# config={}
# check_install()
