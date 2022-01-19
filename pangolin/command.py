#!/usr/bin/env python3
from pangolin import __version__
from . import _program

import os
import sys
import argparse

try:
    import snakemake
except:
    sys.stderr.write(cyan(f'Error: package `{snakemake}` not found, please install snakemake or update pangolin environment.\n'))
    sys.exit(-1)

from tempfile import TemporaryDirectory, TemporaryFile, gettempdir, tempdir
import tempfile
import gzip
import joblib
import select


from pangolin.utils.log_colours import green,cyan
from pangolin.utils import dependency_checks

from pangolin.utils import data_install_checks
from pangolin.utils import update

import pangolin.utils.custom_logger as custom_logger
from pangolin.utils.config import *
from pangolin.utils.initialising import *

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

    a_group = parser.add_argument_group('Analysis modes')
    a_group.add_argument('--accurate', action="store_true",help="Use high accuracy mode for lineage inference. Pipeline: UShER pipeline.")
    a_group.add_argument('--fast', action="store_true",help="Use fast mode for lineage inference. Pipeline: pangoLEARN pipeline.")

    a_group.add_argument('--usher', action="store_true",help="Run UShER pipeline for lineage inference.")
    a_group.add_argument('--pangolearn', action="store_true",help="Run pangoLEARN pipeline for lineage inference.")
    
    a_group.add_argument('--assignment-cache', action="store_true",help="Use cache file from pango-assignment to speed up lineage assignment.", dest="assignment_cache")
    ao_group.add_argument("--no-designation-cache", action='store_true', default=False, help="Developer option - do not use designation hash to assign lineages",dest="no_designation_cache")

    ao_group = parser.add_argument_group('Analysis options')
    ao_group.add_argument('--max-ambig', action="store", default=0.3, type=float,help="Maximum proportion of Ns allowed for pangolin to attempt assignment. Default: 0.3",dest="maxambig")
    ao_group.add_argument('--min-length', action="store", default=25000, type=int,help="Minimum query length allowed for pangolin to attempt assignment. Default: 25000",dest="minlen")

    d_group = parser.add_argument_group('Data options')
    d_group.add_argument("--update", action='store_true', default=False, help="Automatically updates to latest release of pangolin, pangoLEARN and constellations, then exits")
    d_group.add_argument("--update-data", action='store_true',dest="update_data", default=False, help="Automatically updates to latest release of pangoLEARN and constellations, then exits")
    d_group.add_argument('-d', '--datadir', action='store',dest="datadir",help="Data directory minimally containing the pangoLEARN model, header files and UShER tree. Default: Installed pangoLEARN package")
    d_group.add_argument('--usher-tree', action='store', dest='usher_protobuf', help="UShER Mutation Annotated Tree protobuf file to use instead of --usher default from pangoLEARN repository or --datadir")

    m_group = parser.add_argument_group('Misc options')
    m_group.add_argument("--aliases", action='store_true', default=False, help="Print pango-designation alias_key.json and exit")
    m_group.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")
    m_group.add_argument("-pv","--pangoLEARN-version", action='version', version=f"pangoLEARN {pangoLEARN.__version__}",help="show pangoLEARN's version number and exit")
    m_group.add_argument("-dv","--pango-designation-version", action='version', version=f"pango-designation {PANGO_VERSION} used for pangoLEARN and UShER training, alias version {pango_designation.__version__}",help="show pango-designation version number used for training and aliases, then exit")
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

    if args.update:
        update.update({'pangolin': __version__,
                'pangolearn': config[KEY_PANGOLEARN_VERSION],
                'constellations': config[KEY_CONSTELLATIONS_VERSION],
                'scorpio': config[KEY_SCORPIO_VERSION],
                'pango-designation': config[KEY_PANGO_VERSION]
                })

    if args.update_data:
        update.update({'pangolearn': config[KEY_PANGOLEARN_VERSION],
                'constellations': config[KEY_CONSTELLATIONS_VERSION],
                'pango-designation': config[KEY_PANGO_VERSION]}, args.datadir)

    # Parsing analysis mode flags to return one of 'usher', 'pangolearn' or 'assignment_cache'
    config[KEY_ANALYSIS_MODE] = set_up_analysis_mode(args.accurate, args.fast, args.usher, args.pangolearn, args.assignment_cache, config[KEY_ANALYSIS_MODE])
    print(green(f"****\nPangolin running in {config[KEY_ANALYSIS_MODE]} mode.\n****"))

#     config[KEY_DATADIR] = check_datadir(args.datadir)

#     setup_data(config[KEY_DATADIR],config[KEY_PANGOLEARN_VERSION], config)

#     # print(f"Looking in {data_dir} for data files...")

#     if args.aliases:
#         print_alias_file_exit(config[KEY_ALIAS_FILE])

#     if args.all_versions:
#         print(f"pangolin: {__version__}\n"
#               f"pangolearn: {pangoLEARN.__version__}\n"
#               f"constellations: {constellations.__version__}\n"
#               f"scorpio: {scorpio.__version__}\n"
#               f"pango-designation used by pangoLEARN/Usher: {PANGO_VERSION}\n"
#               f"pango-designation aliases: {pango_designation.__version__}")
#         sys.exit(0)

#     # to enable not having to pass a query if running update
#     # by allowing query to accept 0 to many arguments
#     if len(args.query) > 1:
#         print(cyan(f"Error: Too many query (input) fasta files supplied: {args.query}\nPlease supply one only."))
#         parser.print_help()
#         sys.exit(-1)

#     # find the query fasta
#     try:
#         if not os.path.exists(os.path.join(cwd, args.query[0])):
#             if select.select([sys.stdin,],[],[],0.0)[0]:
#                 query = sys.stdin
#             elif not select.select([sys.stdin,],[],[],0.0)[0]:
#                 tried_path = os.path.join(cwd, args.query[0])
#                 if tried_path.endswith("-"):
#                     sys.stderr.write(cyan(
#                         f'Error: cannot find query (input) fasta file using stdin.\n' +
#                                     'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
#                                     ' for detailed instructions.\n'))
#                     sys.exit(-1)
#                 else:
#                     sys.stderr.write(cyan(f'Error: cannot find query (input) fasta file at:') + f'{tried_path}\n' +
#                                     'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
#                                     ' for detailed instructions.\n')
#                     sys.exit(-1)
#         else:
#             query = os.path.join(cwd, args.query[0])
#             print(green(f"The query file is:") + f"{query}")
#     except IndexError:
#         sys.stderr.write(cyan(
#             f'Error: input query fasta could not be detected from a filepath or through stdin.\n' +
#             'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
#             ' for detailed instructions.\n'))
#         sys.exit(-1)

#     # default output dir

#     if args.outdir:
#         outdir = os.path.join(cwd, args.outdir)
#         if not os.path.exists(outdir):
#             try:
#                 os.mkdir(outdir)
#             except:
#                 sys.stderr.write(cyan(f'Error: cannot create directory:') + f"{outdir}")
#                 sys.exit(-1)
#     else:
#         outdir = cwd

#     if args.outfile:
#         outfile = os.path.join(outdir, args.outfile)
#     else:
#         outfile = os.path.join(outdir, "lineage_report.csv")

#     if args.tempdir:
#         to_be_dir = os.path.join(cwd, args.tempdir)
#         if not os.path.exists(to_be_dir):
#             os.mkdir(to_be_dir)
#         temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
#         tempdir = temporary_directory.name
#     else:
#         temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
#         tempdir = temporary_directory.name

#     if args.no_temp:
#         print(green(f"\n--no-temp: ") + f"all intermediate files will be written to {outdir}\n")
#         tempdir = outdir

#     if args.alignment:
#         align_dir = outdir
#         alignment_out = True
#     else:
#         align_dir = tempdir
#         alignment_out = False

#     """
#     QC steps:
#     1) check no empty seqs
#     2) check N content
#     3) write a file that contains just the seqs to run
#     """

#     # do_not_run = []
#     # run = []


#     data_install_checks.check_install(config)
#     snakefile = data_install_checks.get_snakefile(thisdir)

#     set_up_verbosity(config)

#     trained_model = ""
#     header_file = ""
#     designated_hash = ""
#     use_usher = args.usher
#     if args.usher_protobuf:
#         usher_protobuf = os.path.join(cwd, args.usher_protobuf)
#         if not os.path.exists(usher_protobuf):
#             sys.stderr.write('Error: cannot find --usher-tree file at {}\n'.format(usher_protobuf))
#             sys.exit(-1)
#         use_usher = True
#     else:
#         usher_protobuf = ""

#     use_cache = args.use_cache
#     cache = ""
#     if use_cache:
#         try:
#             import pangolin_assignment
#             pangolin_assignment_dir = pangolin_assignment.__path__[0]
#             for r, d, f in os.walk(pangolin_assignment_dir):
#                 for fn in f:
#                     if fn == "pango_assignment.cache.csv.gz" and cache == "":
#                         cache = os.path.join(r, fn)
#             if not os.path.exists(cache):
#                 sys.stderr.write('Error: cannot find cache\n')
#                 sys.exit(-1)
#         except:
#             sys.stderr.write(cyan('Error: please install `pangolin_assignment` with \n') +
#                              "pip install git+https://github.com/cov-lineages/pangolin-assignment.git")
#             sys.exit(-1)

#         try:
#             with gzip.open(cache, 'rt') as f:
#                 line = f.readline()
#         except:
#             with open(cache, 'r') as f:
#                 line = f.readline()
#                 if "git-lfs.github.com" in line:
#                     sys.stderr.write(
#                         'Error: Git LFS file not pulled successfully. Please install git-lfs \nusing conda or an alternative (not pip) then re-install pangolin-assignment \nwith pip install git+https://github.com/cov-lineages/pangolin-assignment.git\n')
#                     sys.exit(-1)

#     for r,d,f in os.walk(data_dir):
#         for fn in f:
#             if fn == "decisionTreeHeaders_v1.joblib":
#                 header_file = os.path.join(r, fn)
#             elif fn == "decisionTree_v1.joblib":
#                 trained_model = os.path.join(r, fn)
#             elif fn =="lineages.hash.csv":
#                 designated_hash = os.path.join(r, fn)
#             elif fn == "lineageTree.pb" and usher_protobuf == "":
#                 usher_protobuf = os.path.join(r, fn)
#     if ((use_usher and (usher_protobuf == "" or designated_hash=="") or
#         (not use_usher and (trained_model=="" or header_file=="" or designated_hash=="")))):
#         print(cyan("""pangoLEARN version should be >= 2021-05-27. \n
# Appropriate data files not found from the installed pangoLEARN repo.
# Please see https://cov-lineages.org/pangolin.html for installation and updating instructions."""))
#         exit(1)

#         print(green("\nData files found:"))
#         if use_usher:
#             print(f"UShER tree:\t{usher_protobuf}")
#             print(f"Designated hash:\t{designated_hash}")
#         else:
#             print(f"Trained model:\t{trained_model}")
#             print(f"Header file:\t{header_file}")
#             print(f"Designated hash:\t{designated_hash}")
#         if use_cache:
#             print(f"Assignment cache:\t{cache}")


#         config["trained_model"] = trained_model
#         config["header_file"] = header_file
#         config["designated_hash"] = designated_hash
#         config["cache"] = cache

#     if use_usher:
#         config["usher_protobuf"] = usher_protobuf

#     if config['verbose']:
#         print(green("\n**** CONFIG ****"))
#         for k in sorted(config):
#             print(green(k), config[k])

#         status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
#                                         workdir=tempdir,config=config, cores=args.threads,lock=False
#                                         )
#     else:
#         logger = custom_logger.Logger()
#         status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=tempdir,
#                                     config=config, cores=args.threads,lock=False,quiet=True,log_handler=config["log_api"]
#                                     )

#     if status: # translate "success" into shell exit code of 0
#        return 0

#     return 1


if __name__ == '__main__':
    main()
