#!/usr/bin/env python3
from . import _program
from pangolin import __version__
from pangolin.utils import data_checks

try:
    import pangolin_data
except:
    data_checks.install_error("pangolin_data", "https://github.com/cov-lineages/pangolin-data.git")

try:
    import scorpio
except:
    data_checks.install_error("scorpio", "https://github.com/cov-lineages/scorpio.git")

try:
    import constellations
except:
    data_checks.install_error("constellations", "https://github.com/cov-lineages/constellations.git")

import os
import sys
import argparse

try:
    import snakemake
except:
    sys.stderr.write(cyan(f'Error: package `{snakemake}` not found, please install snakemake or update pangolin environment.\n'))
    sys.exit(-1)


from pangolin.utils.log_colours import green,cyan
from pangolin.utils import dependency_checks

from pangolin.utils import update


from pangolin.utils.config import *
from pangolin.utils.initialising import *
import pangolin.utils.io_parsing as io

from pangolin.utils.report_collation import generate_final_report,get_voc_list

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()


def main(sysargs = sys.argv[1:]):
    parser = argparse.ArgumentParser(prog = _program,
    description='pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages',
    usage='''pangolin <query> [options]''')

    io_group = parser.add_argument_group('Input-Output options')
    io_group.add_argument('query', nargs="*", help='Query fasta file of sequences to analyse.')
    io_group.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    io_group.add_argument('--outfile', action="store",help="Optional output file name. Default: lineage_report.csv")
    io_group.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    io_group.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    io_group.add_argument('--alignment', action="store_true",help="Output multiple sequence alignment.")
    io_group.add_argument('--alignment-file', action="store",help="Multiple sequence alignment file name.")
    io_group.add_argument('--expanded-lineage', action="store_true", default=False, help="Optional expanded lineage from alias.json in report.")


    a_group = parser.add_argument_group('Analysis options')
    a_group.add_argument('--analysis-mode', action="store",help="Specify which inference engine to use. Options: accurate (UShER), fast (pangoLEARN), pangolearn, usher. Default: UShER inference.")
    
    a_group.add_argument("--skip-designation-cache", action='store_true', default=False, help="Developer option - do not use designation cache to assign lineages.",dest="skip_designation_cache")
    a_group.add_argument("--skip-scorpio", action='store_true', default=False, help="Developer option - do not use scorpio to check VOC/VUI lineage assignments.",dest="skip_scorpio")

    a_group.add_argument('--max-ambig', action="store", default=0.3, type=float,help="Maximum proportion of Ns allowed for pangolin to attempt assignment. Default: 0.3",dest="maxambig")
    a_group.add_argument('--min-length', action="store", default=25000, type=int,help="Minimum query length allowed for pangolin to attempt assignment. Default: 25000",dest="minlen")
    a_group.add_argument('--usher', action='store_true', default=False, help=argparse.SUPPRESS)

    d_group = parser.add_argument_group('Data options')
    d_group.add_argument("--update", action='store_true', default=False, help="Automatically updates to latest release of pangolin, pangolin-data, scorpio and constellations (and pangolin-assignment if it has been installed using --add-assignment-cache), then exits.")
    d_group.add_argument("--update-data", action='store_true',dest="update_data", default=False, help="Automatically updates to latest release of constellations and pangolin-data, including the pangoLEARN model, UShER tree file and alias file (also pangolin-assignment if it has been installed using --add-assignment-cache), then exits.")
    d_group.add_argument('--add-assignment-cache', action='store_true', dest="add_assignment_cache", default=False, help="Install the pangolin-assignment repository for use with --use-assignment-cache.  This makes updates slower and makes pangolin slower for small numbers of input sequences but much faster for large numbers of input sequences.")
    d_group.add_argument('--use-assignment-cache', action='store_true', dest="use_assignment_cache", default=False, help="Use assignment cache from optional pangolin-assignment repository. NOTE: the repository must be installed by --add-assignment-cache before using --use-assignment-cache.")
    d_group.add_argument('-d', '--datadir', action='store',dest="datadir",help="Data directory minimally containing the pangoLEARN model, header files and UShER tree. Default: Installed pangolin-data package.")
    d_group.add_argument('--usher-tree', action='store', dest='usher_protobuf', help="UShER Mutation Annotated Tree protobuf file to use instead of default from pangolin-data repository or --datadir.")
    d_group.add_argument('--assignment-cache', action='store', dest='assignment_cache', help="Cached precomputed assignment file to use instead of default from pangolin-assignment repository.  Does not require installation of pangolin-assignment.")

    m_group = parser.add_argument_group('Misc options')
    m_group.add_argument("--aliases", action='store_true', default=False, help="Print Pango alias_key.json and exit.")
    m_group.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")
    m_group.add_argument("-pv","--pangolin-data-version", action='version', version=f"pangolin-data {pangolin_data.__version__}",help="show version number of pangolin data files (UShER tree and pangoLEARN model files) and exit.")
    m_group.add_argument("--all-versions", action='store_true',dest="all_versions", default=False, help="Print all tool, dependency, and data versions then exit.")
    m_group.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    m_group.add_argument("-t","--threads",action="store",default=1,type=int, help="Number of threads")


    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # Initialise config dict
    config = setup_config_dict(cwd)
    data_checks.check_install(config)
    set_up_verbosity(config)

    if args.usher:
        sys.stderr.write(cyan(f"--usher is a pangolin v3 option and is deprecated in pangolin v4.  UShER is now the default analysis mode.  Use --analysis-mode to explicitly set mode.\n"))

    setup_data(args.datadir,config[KEY_ANALYSIS_MODE], config)

    if args.add_assignment_cache:
        update.install_pangolin_assignment()

    if args.update:
        version_dictionary = {'pangolin': __version__,
                              'pangolin-data': config[KEY_PANGOLIN_DATA_VERSION],
                              'constellations': config[KEY_CONSTELLATIONS_VERSION],
                              'scorpio': config[KEY_SCORPIO_VERSION]}
        update.add_pangolin_assignment_if_installed(version_dictionary)
        update.update(version_dictionary)

    if args.update_data:
        version_dictionary = {'pangolin-data': config[KEY_PANGOLIN_DATA_VERSION],
                              'constellations': config[KEY_CONSTELLATIONS_VERSION]}
        update.add_pangolin_assignment_if_installed(version_dictionary)
        update.update(version_dictionary, args.datadir)

    # install_pangolin_assignment doesn't exit so that --update/--update-data can be given at the
    # same time (or a query file).  If --add-assignment-cache is the only arg, exit without error.
    if args.add_assignment_cache and not args.query:
        sys.exit(0)

    # add flag to config for whether to run scorpio
    if args.skip_scorpio:
        print(green(f"****\nPangolin skipping scorpio steps.\n****"))
        config[KEY_SKIP_SCORPIO] = True
    
    if args.expanded_lineage:
        print(green(f"****\nAdding expanded lineage column to output.\n****"))
        config[KEY_EXPANDED_LINEAGE] = True
        
    # Parsing analysis mode flags to return one of 'usher' or 'pangolearn'
    config[KEY_ANALYSIS_MODE] = set_up_analysis_mode(args.analysis_mode, config[KEY_ANALYSIS_MODE])
    print(green(f"****\nPangolin running in {config[KEY_ANALYSIS_MODE]} mode.\n****"))
    snakefile = get_snakefile(thisdir,config[KEY_ANALYSIS_MODE])

    config[KEY_DESIGNATION_CACHE],config[KEY_ALIAS_FILE] = data_checks.find_designation_cache_and_alias(config[KEY_DATADIR],DESIGNATION_CACHE_FILE,ALIAS_FILE)
    if args.aliases:
        print_alias_file_exit(config[KEY_ALIAS_FILE])

    if args.all_versions:
        print_versions_exit(config)

    # to enable not having to pass a query if running update
    # by allowing query to accept 0 to many arguments
    
#   setup outdir and outfiles
    config[KEY_OUTDIR] = io.set_up_outdir(args.outdir,cwd,config[KEY_OUTDIR])
    config[KEY_OUTFILE] = io.set_up_outfile(args.outfile, config[KEY_OUTFILE],config[KEY_OUTDIR])
    io.set_up_tempdir(args.tempdir,args.no_temp,cwd,config[KEY_OUTDIR], config)
    config[KEY_ALIGNMENT_FILE],config[KEY_ALIGNMENT_OUT] = io.parse_alignment_options(args.alignment, config[KEY_OUTDIR], config[KEY_TEMPDIR],args.alignment_file, config[KEY_ALIGNMENT_FILE])

    config[KEY_QUERY_FASTA] = io.find_query_file(cwd, config[KEY_TEMPDIR], args.query)

    io.quick_check_query_file(cwd, args.query, config[KEY_QUERY_FASTA])

    if config[KEY_ANALYSIS_MODE] == "usher":
        # Find usher protobuf file (and if specified, assignment cache file too)
        data_checks.get_datafiles(config[KEY_DATADIR],usher_files,config)
        if args.usher_protobuf:
            config[KEY_USHER_PB] = data_checks.check_file_arg(args.usher_protobuf, cwd, '--usher-tree')
            print(green(f"Using usher tree file {args.usher_protobuf}"))
        if args.assignment_cache:
            config[KEY_ASSIGNMENT_CACHE] = data_checks.check_file_arg(args.assignment_cache, cwd, '--assignment-cache')
            print(green(f"Using assignment cache file {args.assignment_cache}"))
        elif args.use_assignment_cache:
            config[KEY_ASSIGNMENT_CACHE] = data_checks.get_assignment_cache(USHER_ASSIGNMENT_CACHE_FILE, config)
            print(green("Using pangolin-assignment cache"))
        else:
            config[KEY_ASSIGNMENT_CACHE] = ""

    elif config[KEY_ANALYSIS_MODE] == "pangolearn":
        # find designation cache and the model files
        data_checks.get_datafiles(config[KEY_DATADIR],pangolearn_files,config)
        if args.use_assignment_cache or args.assignment_cache:
            sys.stderr.write(cyan(f"Warning: --use-assignment-cache and --assignment-cache are ignored when --analysis-mode is 'fast' or 'pangolearn'.\n"))

    preprocessing_snakefile = get_snakefile(thisdir,"preprocessing")

    if args.verbose:
        print(green("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(k), config[k])

        status = snakemake.snakemake(preprocessing_snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                        workdir=config[KEY_TEMPDIR],config=config, cores=args.threads,lock=False
                                        )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(preprocessing_snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=config[KEY_TEMPDIR],
                                    config=config, cores=args.threads,lock=False,quiet=True,log_handler=logger.log_handler
                                    )
    if status: # translate "success" into shell exit code of 0
       
        if config[KEY_VERBOSE]:
            print(green("\n**** CONFIG ****"))
            for k in sorted(config):
                print(green(k), config[k])

            status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                            workdir=config[KEY_TEMPDIR],config=config, cores=args.threads,lock=False
                                            )
        else:
            logger = custom_logger.Logger()
            status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=config[KEY_TEMPDIR],
                                        config=config, cores=args.threads,lock=False,quiet=True,log_handler=logger.log_handler
                                        )
        
       
        if status: 
            
            ## Collate the report here

            preprocessing_csv = os.path.join(config[KEY_TEMPDIR],"preprocessing.csv")
            inference_csv = os.path.join(config[KEY_TEMPDIR],"inference_report.csv")
            cached_csv = os.path.join(config[KEY_TEMPDIR],"cache_assigned.csv")
            constellation_list = get_voc_list(os.path.join(config[KEY_TEMPDIR], "get_constellations.txt"), config[KEY_ALIAS_FILE])

            generate_final_report(preprocessing_csv, inference_csv, cached_csv, config[KEY_ALIAS_FILE], constellation_list, config[KEY_PANGOLIN_DATA_VERSION],config[KEY_ANALYSIS_MODE], args.skip_designation_cache, config[KEY_OUTFILE],config)

            print(green(f"****\nOutput file written to: ") + config[KEY_OUTFILE])

            if config[KEY_ALIGNMENT_OUT]:
                print(green(f"****\nOutput alignment written to: ") + config[KEY_ALIGNMENT_FILE])


            return 0

        return 1
    return 1

if __name__ == '__main__':
    main()
