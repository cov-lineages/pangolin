from pangolin.utils.log_colours import green,cyan
from pangolin.utils.config import *
from pangolin.utils.error_messages import *

try:
    import pangoLEARN
except:
    install_error("pangoLEARN", "https://github.com/cov-lineages/pangoLEARN.git")

try:
    from pangoLEARN import PANGO_VERSION
except:
    sys.stderr.write(cyan('Error: please update to pangoLEARN version >= 2021-05-27\n'))
    sys.exit(-1)

try:
    import scorpio
except:
    install_error("scorpio", "https://github.com/cov-lineages/scorpio.git")

try:
    import constellations
except:
    install_error("constellations", "https://github.com/cov-lineages/constellations.git")

try:
    import pango_designation
except:
    install_error("pango_designation", "https://github.com/cov-lineages/pango-designation.git")


def setup_config_dict(cwd):
    default_dict = {
            KEY_ANALYSIS_MODE:"accurate", #options: accurate, fast, usher, pangolearn, cache-assign
            
            KEY_SKIP_DESIGNATION_HASH: False,
            KEY_USE_CACHE: False,

            KEY_QUERY_FASTA:None,

            KEY_OUTDIR:cwd,
            KEY_OUTFILE:"lineage_report.csv",

            KEY_ALIGNDIR: cwd,
            KEY_ALIGNMENT_OUT: False,

            KEY_ALIGNMENT_FILE:"alignment.fasta",

            KEY_TEMPDIR:None,
            KEY_NO_TEMP:False,
            
            KEY_DATADIR:None,

            KEY_TRIM_START:265, # where to pad to using datafunk
            KEY_TRIM_END:29674, # where to pad after using datafunk
            
            KEY_ALIAS_FILE: None,

            KEY_CONSTELLATION_FILES: constellation_files,
            
            KEY_PANGOLEARN_VERSION: pangoLEARN.__version__,
            KEY_PANGOLIN_VERSION: __version__,
            KEY_PANGO_VERSION: PANGO_VERSION,
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
    
    if usher_arg:
        accurate_arg = True
    if pangolearn_arg:
        fast_arg = True

    analysis_options = {"accurate":accurate_arg,
                        "fast":fast_arg,
                        "cache":assignment_cache_arg}
    option_count = 0
    
    for i in analysis_options:
        if analysis_options[i]==True:
            option_count +=1
            analysis_mode = i

    if option_count > 1:
        sys.stderr.write(cyan(f"Incompatible options specified: please select one of `--fast`,`--accurate`,`--pangolearn`,`--usher` or `--assignment-cache`\n"))
        sys.exit(-1)
    
    return analysis_mode
    

def check_datadir(datadir_arg):

    # find the data
    if datadir_arg:
        # this needs to be an absolute path when we pass it to scorpio
        datadir = os.path.abspath(args.datadir)
        if not os.path.exists(config[KEY_DATADIR]):
            sys.stderr.write(cyan(f"Cannot find data directory specified: {config[KEY_DATADIR]}\n"))
            sys.exit(-1)
    return datadir


def setup_data(datadir,pangolearn_version, config)
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
        for r,d,f in datadir:
            for fn in f:
                if r.endswith('pangoLEARN') and fn == "__init__.py":
                    # print("Found __init__.py")
                    version = version_from_init(os.path.join(r, fn))
                    if version > pangolearn_version:
                        # only use this for pangoLEARN if the version is > than what we already have
                        pangolearn_version = version
                        use_datadir = True
                        
    if use_datadir == False:
        # we haven't got a viable datadir from searching args.datadir
        pangoLEARN_dir = pangoLEARN.__path__[0]
        datadir = os.path.join(pangoLEARN_dir,"data")

    if not alias_file:
        sys.stderr.write(cyan('Error: Could not find alias file'))
        install_error("pango-designation", "https://github.com/cov-lineages/pango-designation.git")

    config[KEY_PANGOLEARN_VERSION] = pangolearn_version
    config[KEY_PANGO_VERSION] = pango_designation_version
    config[KEY_CONSTELLATIONS_VERSION] = constellations_version
    config[KEY_ALIAS_FILE] = alias_file
    config[KEY_DATADIR] = datadir

def print_alias_file_exit(alias_file):
    with open(alias_file, 'r') as handle:
        for line in handle:
            print(line.rstrip())
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

        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "