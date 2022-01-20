#!/usr/bin/env python3
from . import _program
from pangolin import __version__

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
    import pango_designation
except:
    install_error("pango-designation", "https://github.com/cov-lineages/pango-designation.git")

try:
    import scorpio
except:
    install_error("scorpio", "https://github.com/cov-lineages/scorpio.git")

try:
    import constellations
except:
    install_error("constellations", "https://github.com/cov-lineages/constellations.git")

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

from pangolin.utils import data_checks
from pangolin.utils import update


from pangolin.utils.config import *
from pangolin.utils.initialising import *
import pangolin.utils.io_parsing as io

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
    a_group.add_argument("--skip-designation-cache", action='store_true', default=False, help="Developer option - do not use designation hash to assign lineages.",dest="skip_designation_cache")

    ao_group = parser.add_argument_group('Analysis options')
    ao_group.add_argument('--max-ambig', action="store", default=0.3, type=float,help="Maximum proportion of Ns allowed for pangolin to attempt assignment. Default: 0.3",dest="maxambig")
    ao_group.add_argument('--min-length', action="store", default=25000, type=int,help="Minimum query length allowed for pangolin to attempt assignment. Default: 25000",dest="minlen")

    d_group = parser.add_argument_group('Data options')
    d_group.add_argument("--update", action='store_true', default=False, help="Automatically updates to latest release of pangolin, pangoLEARN and constellations, then exits.")
    d_group.add_argument("--update-data", action='store_true',dest="update_data", default=False, help="Automatically updates to latest release of pangoLEARN and constellations, then exits.")
    d_group.add_argument('-d', '--datadir', action='store',dest="datadir",help="Data directory minimally containing the pangoLEARN model, header files and UShER tree. Default: Installed pangoLEARN package.")
    d_group.add_argument('--usher-tree', action='store', dest='usher_protobuf', help="UShER Mutation Annotated Tree protobuf file to use instead of --usher default from pangoLEARN repository or --datadir.")

    m_group = parser.add_argument_group('Misc options')
    m_group.add_argument("--aliases", action='store_true', default=False, help="Print pango-designation alias_key.json and exit.")
    m_group.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")
    m_group.add_argument("-pv","--pangoLEARN-version", action='version', version=f"pangoLEARN {pangoLEARN.__version__}",help="show pangoLEARN's version number and exit.")
    m_group.add_argument("-dv","--pango-designation-version", action='version', version=f"pango-designation {PANGO_VERSION} used for pangoLEARN and UShER training, alias version {pango_designation.__version__}",help="show pango-designation version number used for training and aliases, then exit.")
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
    snakefile = get_snakefile(thisdir,config[KEY_ANALYSIS_MODE])

    setup_data(args.datadir,config[KEY_ANALYSIS_MODE], config)

    if args.aliases:
        print_alias_file_exit(config[KEY_ALIAS_FILE])

    if args.all_versions:
        print_versions_exit(config)

    # to enable not having to pass a query if running update
    # by allowing query to accept 0 to many arguments
    config[KEY_QUERY_FASTA] = io.find_query_file(cwd, args.query)

    io.quick_check_query_file(config[KEY_QUERY_FASTA])
#   setup outdir and outfiles
    config[KEY_OUTDIR] = io.set_up_outdir(args.outdir,cwd,config[KEY_OUTDIR])
    config[KEY_OUTFILE] = io.set_up_outfile(args.outfile, config[KEY_OUTFILE],config[KEY_OUTDIR])
    config[KEY_TEMPDIR] = io.set_up_tempdir(args.tempdir,args.no_temp,cwd,config[KEY_OUTDIR])
    config[KEY_ALIGNMENT_FILE],config[KEY_ALIGNMENT_OUT] = io.parse_alignment_options(args.alignment, config[KEY_OUTDIR], config[KEY_TEMPDIR],args.alignment_file, config[KEY_ALIGNMENT_FILE])
    
    config[KEY_DESIGNATION_CACHE] = data_checks.find_designation_cache(config[KEY_DATADIR],designation_cache_file,args.skip_designation_cache)

    if config[KEY_ANALYSIS_MODE] == "usher":
        # needed data is usher protobuf file
        config[KEY_USHER_PB] = data_checks.get_usher_protobuf_arg(args.usher_protobuf,cwd)
        data_checks.get_datafiles(config[KEY_DATADIR],usher_files,config)

    elif config[KEY_ANALYSIS_MODE] == "pangolearn":
        # find designation cache and the model files
        data_checks.get_datafiles(config[KEY_DATADIR],pangolearn_files,config)

    elif config[KEY_ANALYSIS_MODE] == "assignment_cache":
        # look for the assignment cache, and also the ??? files (usher or pangolearn?)
        config[KEY_ASSIGNMENT_CACHE] = data_checks.get_cache()

#  """
#     QC steps:
#     1) check no empty seqs
#     2) check N content
#     3) write a file that contains just the seqs to run
#     """
#         # do_not_run = []
#         # run = []
        
#         print(green("** Running sequence QC **"))

#         if os.path.exists(os.path.join(cwd, args.query[0])):
#             file_ending = query.split(".")[-1]
#             if file_ending in ["gz","gzip","tgz"]:
#                 query = gzip.open(query, 'rt')
#             elif file_ending in ["xz","lzma"]:
#                 query = lzma.open(query, 'rt')
                
#         post_qc_query = os.path.join(tempdir, 'query.post_qc.fasta')
#         fw_pass = open(post_qc_query,"w")
#         qc_fail = os.path.join(tempdir,'query.failed_qc.fasta')
#         fw_fail = open(qc_fail,"w")

#         total_input = 0
#         total_pass = 0
        
#         try:
#             for record in SeqIO.parse(query, "fasta"):
#                 total_input +=1
#                 record.description = record.description.replace(' ', '_').replace(",","_")
#                 record.id = record.description
#                 if "," in record.id:
#                     record.id=record.id.replace(",","_")

#                 if len(record) <args.minlen:
#                     record.description = record.description + f" fail=seq_len:{len(record)}"
#                     fw_fail.write(f">{record.description}\n{record.seq}\n")
#                 else:
#                     num_N = str(record.seq).upper().count("N")
#                     prop_N = round((num_N)/len(record.seq), 2)
#                     if prop_N > args.maxambig:
#                         record.description = record.description + f" fail=N_content:{prop_N}"
#                         fw_fail.write(f">{record.description}\n{record.seq}\n")
#                     else:
#                         total_pass +=1
#                         seq = str(record.seq).replace("-","")
#                         fw_pass.write(f">{record.description}\n{seq}\n")
#         except UnicodeDecodeError:
#             sys.stderr.write(cyan(
#                 f'Error: the input query fasta could not be parsed.\n' +
#                 'Double check your query fasta and that compressed stdin was not passed.\n' +
#                 'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
#                 ' for detailed instructions.\n'))
#             sys.exit(-1)

#         print(green("Number of sequences detected: ") + f"{total_input}")
#         print(green("Total passing QC: ") + f"{total_pass}")
#         fw_fail.close()
#         fw_pass.close()

#         if total_pass == 0:
#             with open(outfile, "w") as fw:
#                 fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
#                 for record in SeqIO.parse(os.path.join(tempdir,'query.failed_qc.fasta'), "fasta"):
#                     desc = record.description.split(" ")
#                     reason = ""
#                     for item in desc:
#                         if item.startswith("fail="):
#                             reason = item.split("=")[1]
#                     fw.write(f"{record.id},None,,,,,,PANGO-{PANGO_VERSION},{__version__},{pangoLEARN.__version__},{PANGO_VERSION},fail,{reason}\n")
#             print(cyan(f'Note: no query sequences have passed the qc\n'))
#             sys.exit(0)

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
