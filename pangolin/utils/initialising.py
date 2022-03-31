import os
import sys
import itertools
from distutils.version import LooseVersion

import pangolin.utils.custom_logger as custom_logger
from pangolin.utils.log_colours import green,cyan
from pangolin.utils.config import *
from pangolin.utils.data_checks import *
from pangolin import __version__

import pangolin_data
import scorpio
import constellations

def setup_config_dict(cwd):
    default_dict = {
            KEY_ANALYSIS_MODE:"usher", #options: accurate, fast, usher, pangolearn
            
            KEY_DESIGNATION_CACHE: "",

            KEY_QUERY_FASTA:None,

            KEY_OUTDIR:cwd,
            KEY_OUTFILE:"lineage_report.csv",

            KEY_ALIGNDIR: cwd,
            KEY_ALIGNMENT_FILE:"alignment.fasta",
            KEY_ALIGNMENT_OUT: False,

            KEY_TEMPDIR:None,
            KEY_NO_TEMP:False,
            
            KEY_DATADIR:None,

            KEY_MINLEN: 25000,
            KEY_MAXAMBIG: 0.3,
            KEY_TRIM_START:265, # where to pad to using datafunk
            KEY_TRIM_END:29674, # where to pad after using datafunk
            
            KEY_ALIAS_FILE: None,

            KEY_SKIP_SCORPIO: False,

            KEY_EXPANDED_LINEAGE: False,

            KEY_CONSTELLATION_FILES: [],
            
            KEY_PANGOLIN_VERSION: __version__,
            KEY_PANGOLIN_DATA_VERSION: pangolin_data.__version__,
            KEY_SCORPIO_VERSION: scorpio.__version__,
            KEY_CONSTELLATIONS_VERSION: constellations.__version__,

            KEY_VERBOSE: False,
            KEY_LOG_API: "",
            KEY_THREADS: 1
            }
    return default_dict

def set_up_analysis_mode(analysis_arg, default_mode):
    """
    the logic here 
    - takes the default mode set in the config dict (accurate)
    - it equates the usher arg to accurate arg and pangolearn to fast
    - checks if incompatible flags were used (only one of accurate, fast or cache)
    - overwrites default if any other analysis mode flagged
    - returns new analysis mode
    """
    
    analysis_mode = default_mode
    if analysis_arg:
        if not analysis_arg in ["usher","pangolearn","fast","accurate"]:
            sys.stderr.write(cyan(f"Invalid `--analysis-mode` option specified: please select one of `fast`,`accurate`,`pangolearn` or`usher`\n"))
            sys.exit(-1)

        if analysis_arg in ['pangolearn','fast']:
            analysis_mode = "pangolearn"
        elif analysis_arg in ['usher','accurate']:
            analysis_mode = "usher"

    return analysis_mode
    
def get_snakefile(thisdir,analysis_mode):
    # in this case now, the snakefile used should be the name of the analysis mode (i.e. pangolearn, usher or preprocessing)
    snakefile = os.path.join(thisdir, 'scripts',f'{analysis_mode}.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write(cyan(f'Error: cannot find Snakefile at {snakefile}. Check installation\n'))
        sys.exit(-1)
    return snakefile

def check_datadir(datadir_arg):
    datadir = None
    # find the data
    if datadir_arg:
        # this needs to be an absolute path when we pass it to scorpio
        datadir = os.path.abspath(datadir_arg)
        if not os.path.exists(datadir):
            sys.stderr.write(cyan(f"Cannot find data directory specified: {datadir}\n"))
            sys.exit(-1)
    return datadir

def version_from_init(init_file):
    version=None
    with open(init_file, "r") as fr:
        for l in fr:
            if l.startswith("__version__"):
                l = l.rstrip("\n")
                version = l.split('=')[1]
                version = version.replace('"',"").replace(" ","")
                break
    return version

def setup_data(datadir_arg,analysis_mode, config):

    datadir = check_datadir(datadir_arg)

    pangolin_data_dir = pangolin_data.__path__[0]
    constellations_dir = constellations.__path__[0]
    constellation_files = []

    data_locations = [os.walk(constellations_dir)]

    if datadir:
        data_locations.append(os.walk(datadir))

    # the logic of this is to search the "built-in" constellations
    # path first and then if as custom datadir is passed, follow up with those, so that
    # any files found in the datadir supercede the "built-in" modules. The assumption
    # here is that the datadir contains newer (user updated) data
    for r, _, f in itertools.chain.from_iterable(data_locations):
        if r.endswith('/constellations') or r.endswith('/constellations/definitions'):
            constellation_files = []  # only collect the constellations from the last directory found
        for fn in f:
            if r.endswith('/constellations') and fn == '__init__.py':
                constellations_version = version_from_init(os.path.join(r, fn))
            elif (r.endswith('/constellations') or r.endswith('/constellations/definitions')) and fn.endswith('.json'):
                constellation_files.append(os.path.join(r, fn))

    pangolin_data_version = pangolin_data.__version__
    use_datadir = False
    datadir_too_old = False
    if datadir:
        version = "Unknown"
        for r,d,f in os.walk(datadir):
            for fn in f:
                # pangolin-data/__init__.py not constellations/__init__.py:
                if r.endswith('data') and fn == "__init__.py":
                    # print("Found " + os.path.join(r, fn))
                    version = version_from_init(os.path.join(r, fn))
                    if not version:
                        continue
                    
                    if LooseVersion(version) >= LooseVersion(pangolin_data.__version__):
                        # only use this if the version is >= than what we already have
                        pangolin_data_version = version
                        use_datadir = True
                    else:
                        datadir_too_old = True
                        sys.stderr.write(cyan(f"Warning: Ignoring specified datadir {datadir} - it contains pangoLEARN model files older ({version}) than those installed ({pangolin_data.__version__})\n"))

    if use_datadir == False:
        # we haven't got a viable datadir from searching args.datadir
        if datadir and not datadir_too_old:
            sys.stderr.write(cyan(
                f"Warning: Ignoring specified datadir {datadir} - could not find __init__.py file to check versions \n"))

        pangolin_data_dir = pangolin_data.__path__[0]
        datadir = os.path.join(pangolin_data_dir,"data")

    config[KEY_PANGOLIN_DATA_VERSION] = pangolin_data_version
    config[KEY_CONSTELLATIONS_VERSION] = constellations_version
    config[KEY_DATADIR] = datadir
    config[KEY_CONSTELLATION_FILES] = constellation_files

    
def print_alias_file_exit(alias_file):
    with open(alias_file, 'r') as handle:
        for line in handle:
            print(line.rstrip())
    
    sys.exit(0)

def print_versions_exit(config):
    print(f"pangolin: {config[KEY_PANGOLIN_VERSION]}\n"
            f"pangolin-data: {config[KEY_PANGOLIN_DATA_VERSION]}\n"
            f"constellations: {config[KEY_CONSTELLATIONS_VERSION]}\n"
            f"scorpio: {config[KEY_SCORPIO_VERSION]}")
    # Report pangolin_assignment version if it is installed, otherwise ignore
    try:
        import pangolin_assignment
        print(f"pangolin-assignment: {pangolin_assignment.__version__}")
    except:
        pass
    sys.exit(0)

def set_up_verbosity(config):
    if config[KEY_VERBOSE]:
        config["quiet"] = False
        config[KEY_LOG_API] = ""
        config["log_string"] = ""
    else:
        config["quiet"] = True
        logger = custom_logger.Logger()
        config[KEY_LOG_API] = logger.log_handler
