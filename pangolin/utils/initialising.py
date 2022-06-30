import os
import sys
import itertools
import re
import subprocess
from distutils.version import LooseVersion
from Bio import SeqIO

import pangolin.utils.custom_logger as custom_logger
from pangolin.utils.log_colours import green,cyan
from pangolin.utils.config import *
from pangolin.utils.data_checks import *
from pangolin import __version__

import pangolin_data
class PangolinAssignmentWrapper():
    __version__ = None
    __path__ = [None]
try:
    import pangolin_assignment
except ImportError:
    # if we can't import the module, leave the variables we replace it with a mock with suitable attributes
    pangolin_assignment = PangolinAssignmentWrapper()
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

            KEY_MAXAMBIG: 0.3,
            KEY_TRIM_START:265, # where to pad to using datafunk
            KEY_TRIM_END:29674, # where to pad after using datafunk
            
            KEY_ALIAS_FILE: None,

            KEY_SKIP_SCORPIO: False,

            KEY_EXPANDED_LINEAGE: False,

            KEY_CONSTELLATION_FILES: [],

            KEY_INPUT_COMPRESSION_TYPE: "plaintext",
            
            KEY_PANGOLIN_VERSION: __version__,
            KEY_PANGOLIN_DATA_VERSION: pangolin_data.__version__,
            KEY_SCORPIO_VERSION: scorpio.__version__,
            KEY_CONSTELLATIONS_VERSION: constellations.__version__,
            KEY_PANGOLIN_ASSIGNMENT_VERSION: pangolin_assignment.__version__,
            KEY_PANGOLIN_ASSIGNMENT_PATH: pangolin_assignment.__path__[0],

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
        if not analysis_arg in ["usher","pangolearn","fast","accurate","scorpio"]:
            sys.stderr.write(cyan(f"Invalid `--analysis-mode` option specified: please select one of `fast`,`accurate`,`pangolearn`, `usher` or `scorpio`\n"))
            sys.exit(-1)

        if analysis_arg in ['pangolearn','fast']:
            analysis_mode = "pangolearn"
        elif analysis_arg in ['usher','accurate']:
            analysis_mode = "usher"
        elif analysis_arg == "scorpio":
            analysis_mode = "scorpio"

    return analysis_mode
    
def get_snakefile(thisdir,analysis_mode):
    snakefile = ""
    if analysis_mode != "scorpio":
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

def setup_data(datadir_arg, analysis_mode, config, use_old_data):
    datadir = check_datadir(datadir_arg)

    config[KEY_PANGOLIN_DATA_VERSION] = pangolin_data.__version__
    config[KEY_DATADIR] = pangolin_data.__path__[0]
    config[KEY_CONSTELLATIONS_VERSION] = constellations.__version__
    config[KEY_CONSTELLATION_FILES] = get_constellation_files(constellations.__path__[0])
    config[KEY_PANGOLIN_ASSIGNMENT_VERSION] = pangolin_assignment.__version__
    config[KEY_PANGOLIN_ASSIGNMENT_PATH] = pangolin_assignment.__path__[0]
   
    if datadir:
        for module_name in ('constellations', 'pangolin_data', 'pangolin_assignment'):
            for r, _, f in os.walk(datadir):
                for fn in f:
                    if r.endswith('/' + module_name) and fn == '__init__.py':
                        version = version_from_init(os.path.join(r, fn))
                        # module_name has been imported so exists in global namespace
                        current_version = getattr(globals()[module_name], '__version__', '0')
                        if use_old_data or current_version is None or LooseVersion(version) >= LooseVersion(current_version):
                            if module_name == "pangolin_data":
                                config[KEY_PANGOLIN_DATA_VERSION] = version
                                config[KEY_DATADIR] = os.path.join(datadir, r)
                            elif module_name == "pangolin_assignment":
                                config[KEY_PANGOLIN_ASSIGNMENT_VERSION] = version
                                config[KEY_PANGOLIN_ASSIGNMENT_PATH] = os.path.join(datadir, r)
                            elif module_name == "constellations":
                                config[KEY_CONSTELLATIONS_VERSION] = version
                                config[KEY_CONSTELLATION_FILES] = get_constellation_files(r)
                        else:
                            sys.stderr.write(cyan(f"Warning: Ignoring {module_name} in specified datadir {datadir} - it contains {module_name} with older ({version}) than those installed ({current_version})\n"))

def parse_qc_thresholds(maxambig, minlen, reference_fasta, config):
    
    if maxambig:
        maxambig = float(maxambig)
        if maxambig <=1 and maxambig >= 0:
            config[KEY_MAXAMBIG] = maxambig
        else:
            sys.stderr.write(cyan(f'Error: `--max-ambiguity` should be a float between 0 and 1.\n'))
            sys.exit(-1)

    if minlen:
        minlen = float(minlen)
        reflen = 0
        for record in SeqIO.parse(reference_fasta,"fasta"):
            reflen = len(record)
        
        if minlen>reflen:
            sys.stderr.write(cyan(f'Error: `--min-length` should be less than the length of the reference: {reflen}.\n'))
            sys.exit(-1)
        else:
            new_maxambig = round(1-(minlen/reflen), 3)
            print(f"Converting minimum length of {minlen} to maximum ambiguity of {new_maxambig}.")
            if new_maxambig < config[KEY_MAXAMBIG]:
                config[KEY_MAXAMBIG] = new_maxambig
        
    print(green(f"Maximum ambiguity allowed is {config[KEY_MAXAMBIG]}.\n****"))

def print_ram_warning(analysis_mode):
    if analysis_mode == "pangolearn":
        print(cyan("Warning: pangoLEARN mode may use a significant amount of RAM, be aware that it will not suit every system."))

def print_alias_file_exit(alias_file):
    with open(alias_file, 'r') as handle:
        for line in handle:
            print(line.rstrip())
    
    sys.exit(0)

def get_version(programs):
    for program in programs:
        cmd = [f"{program} --version"]
        output = subprocess.run(cmd, shell=True, check=True,
                                stdout=subprocess.PIPE, encoding='utf-8')
        version = output.stdout.strip().split()[-1].strip('()v')
        print(f"{program.split()[0]} {version}")

def print_conda_version(pkg_list):
    for pkg in pkg_list:
        try:
            result = subprocess.run(['conda', 'list', pkg],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
            
        except subprocess.CalledProcessError as e:
            stderr = e.stderr.decode('utf-8')
            sys.stderr.write(cyan(f"Error: {e}:\n{stderr}\n"))
            sys.exit(-1)
        output = result.stdout.decode('utf-8')
        m = re.search(f'\n{pkg} +([0-9.]+) .*\sbioconda$', output)
        if m:
            version = m.group(1)
            print(f"{pkg}: {version}")
        else:
            sys.stderr.write(cyan(f"version not found in output of 'conda list {pkg}':\n{output}\n"))

def print_versions_exit(config):
    print(f"pangolin: {config[KEY_PANGOLIN_VERSION]}\n"
            f"pangolin-data: {config[KEY_PANGOLIN_DATA_VERSION]}\n"
            f"constellations: {config[KEY_CONSTELLATIONS_VERSION]}\n"
            f"scorpio: {config[KEY_SCORPIO_VERSION]}")
    # Report pangolin_assignment version if it is installed, otherwise ignore
    if config[KEY_PANGOLIN_ASSIGNMENT_VERSION] is not None:
        print(f"pangolin-assignment: {config[KEY_PANGOLIN_ASSIGNMENT_VERSION]}")
    # Print versions of other important tools used by pangolin
    get_version(['usher', 'gofasta', 'minimap2'])
    # print_conda_version(['usher', 'ucsc-fatovcf', 'gofasta', 'minimap2'])
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
