import os
import sys
import pangolin.utils.custom_logger as custom_logger
import itertools

from pangolin.utils.log_colours import green,cyan
from pangolin.utils.config import *
from pangolin.utils.error_messages import *
from pangolin import __version__

import pangoLEARN
from pangoLEARN import PANGO_VERSION
import scorpio
import constellations
import pango_designation

def setup_config_dict(cwd):
    default_dict = {
            KEY_ANALYSIS_MODE:"usher", #options: accurate, fast, usher, pangolearn, assignment_cache
            
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

            KEY_TRIM_START:265, # where to pad to using datafunk
            KEY_TRIM_END:29674, # where to pad after using datafunk
            
            KEY_ALIAS_FILE: None,

            KEY_CONSTELLATION_FILES: None,
            
            KEY_PANGOLEARN_VERSION: pangoLEARN.__version__,
            KEY_PANGOLIN_VERSION: __version__,
            KEY_PANGO_VERSION: PANGO_VERSION,
            KEY_PANGO_DESIGNATION_VERSION: pango_designation.__version__,
            KEY_SCORPIO_VERSION: scorpio.__version__,
            KEY_CONSTELLATIONS_VERSION: constellations.__version__,

            KEY_VERBOSE: False,
            KEY_THREADS: 1
            }
    return default_dict

def set_up_analysis_mode(accurate_arg, fast_arg, usher_arg, pangolearn_arg, assignment_cache_arg,default_mode):
    """
    the logic here 
    - takes the default mode set in the config dict (accurate)
    - it equates the usher arg to accurate arg and pangolearn to fast
    - checks if incompatible flags were used (only one of accurate, fast or cache)
    - overwrites default if any other analysis mode flagged
    - returns new analysis mode
    """
    
    analysis_mode = default_mode
    
    if accurate_arg:
        usher_arg = True
    if fast_arg:
        pangolearn_arg = True

    analysis_options = {"usher":usher_arg,
                        "pangolearn":pangolearn_arg,
                        "assignment_cache":assignment_cache_arg}
    option_count = 0
    
    for i in analysis_options:
        if analysis_options[i]==True:
            option_count +=1
            analysis_mode = i

    if option_count > 1:
        sys.stderr.write(cyan(f"Incompatible options specified: please select one of `--fast`,`--accurate`,`--pangolearn`,`--usher` or `--assignment-cache`\n"))
        sys.exit(-1)
    
    return analysis_mode
    
def get_snakefile(thisdir,analysis_mode):
    # in this case now, the snakefile used should be the name of the analysis mode (i.e. pangolearn, usher or cache)
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
    with open(init_file, "r") as fr:
        for l in fr:
            if l.startswith("__version__"):
                l = l.rstrip("\n")
                version = l.split('=')[1]
                version = version.replace('"',"").replace(" ","")
                break
    return version

def pango_version_from_init(init_file):
    with open(init_file, "r") as fr:
        for l in fr:
            if l.startswith("PANGO_VERSION"):
                l = l.rstrip("\n")
                version = l.split('=')[1]
                version = version.replace('"',"").replace(" ","")
                break
    return version


def setup_data(datadir_arg,analysis_mode, config):

    datadir = check_datadir(datadir_arg)

    pango_designation_dir = pango_designation.__path__[0]
    constellations_dir = constellations.__path__[0]
    constellation_files = []

    data_locations = [os.walk(pango_designation_dir), os.walk(constellations_dir)]

    if datadir:
        data_locations.append(os.walk(datadir))

    # the logic of this is to search the "built-in" pango_designation and constellations
    # paths first and then if as custom datadir is passed, follow up with those, so that
    # any files found in the datadir supercede the "built-in" modules. The assumption
    # here is that the datadir contains newer (user updated) data
    for r, _, f in itertools.chain.from_iterable(data_locations):
        if r.endswith('/constellations') or r.endswith('/constellations/definitions'):
            constellation_files = []  # only collect the constellations from the last directory found
        for fn in f:
            if r.endswith('/pango_designation') and fn == "alias_key.json":
                alias_file = os.path.join(r, fn)
                # the __init__.py file for pango_designation is on the same level as alias_key.json
                pango_designation_version = version_from_init(os.path.join(r, '__init__.py'))
            elif r.endswith('/constellations') and fn == '__init__.py':
                constellations_version = version_from_init(os.path.join(r, fn))
            elif (r.endswith('/constellations') or r.endswith('/constellations/definitions')) and fn.endswith('.json'):
                constellation_files.append(os.path.join(r, fn))

    use_datadir = False
    if datadir:
        version = "Unknown"
        for r,d,f in os.walk(datadir):
            for fn in f:
                if r.endswith('pangoLEARN') and fn == "__init__.py":
                    # print("Found __init__.py")
                    version = version_from_init(os.path.join(r, fn))
                    
                    if version > pangoLEARN.__version__:
                        # only use this for pangoLEARN if the version is > than what we already have
                        pangoLEARN.__version__ = version
                        use_datadir = True
                        pango_version = pango_version_from_init(init_file)
    
    if not alias_file:
        sys.stderr.write(cyan('Error: Could not find alias file'))
        install_error("pango-designation", "https://github.com/cov-lineages/pango-designation.git")

    if use_datadir == False:
        # we haven't got a viable datadir from searching args.datadir
        pangoLEARN_dir = pangoLEARN.__path__[0]
        datadir = os.path.join(pangoLEARN_dir,"data")
    else:
        config[PANGO_VERSION] = pango_version

    config[KEY_PANGOLEARN_VERSION] = pangoLEARN.__version__
    config[KEY_PANGO_DESIGNATION_VERSION] = pango_designation_version
    
    config[KEY_CONSTELLATIONS_VERSION] = constellations_version
    config[KEY_ALIAS_FILE] = alias_file
    config[KEY_DATADIR] = datadir

    
def print_alias_file_exit(alias_file):
    with open(alias_file, 'r') as handle:
        for line in handle:
            print(line.rstrip())
        sys.exit(0)

def print_versions_exit(config):
    print(f"pangolin: {config[KEY_PANGOLIN_VERSION]}\n"
            f"pangolearn: {config[KEY_PANGOLEARN_VERSION]}\n"
            f"constellations: {config[KEY_CONSTELLATIONS_VERSION]}\n"
            f"scorpio: {config[KEY_SCORPIO_VERSION]}\n"
            f"pango-designation used by pangoLEARN/Usher: {config[KEY_PANGO_VERSION]}\n"
            f"pango-designation aliases: {config[KEY_PANGO_DESIGNATION_VERSION]}")
    sys.exit(0)

def set_up_verbosity(config):
    if config[KEY_VERBOSE]:
        config["quiet"] = False
        config["log_api"] = ""
        config["log_string"] = ""
    else:
        config["quiet"] = True
        logger = custom_logger.Logger()
        config["log_api"] = logger.log_handler
