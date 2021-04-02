#!/usr/bin/env python3
from pangolin import __version__
import argparse
import os.path
import snakemake
import sys
from urllib import request
from distutils.version import LooseVersion
import subprocess
import json
from tempfile import gettempdir
import tempfile
import pprint
import json
import os
import joblib
import pangoLEARN

import custom_logger as custom_logger
import log_handler_handle as lh
import pangofunks as pfunk

import pkg_resources
from Bio import SeqIO


from . import _program


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program,
    description='pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages',
    usage='''pangolin <query> [options]''')

    parser.add_argument('query', nargs="*", help='Query fasta file of sequences to analyse.')
    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('--outfile', action="store",help="Optional output file name. Default: lineage_report.csv")
    parser.add_argument('--alignment', action="store_true",help="Optional alignment output.")
    parser.add_argument('-d', '--datadir', action='store',dest="datadir",help="Data directory minimally containing a fasta alignment and guide tree")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    parser.add_argument('--decompress-model',action="store_true",dest="decompress",help="Permanently decompress the model file to save time running pangolin.")
    parser.add_argument('--max-ambig', action="store", default=0.5, type=float,help="Maximum proportion of Ns allowed for pangolin to attempt assignment. Default: 0.5",dest="maxambig")
    parser.add_argument('--min-length', action="store", default=10000, type=int,help="Minimum query length allowed for pangolin to attempt assignment. Default: 10000",dest="minlen")
    parser.add_argument('--panGUIlin', action='store_true',help="Run web-app version of pangolin",dest="panGUIlin")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument("-t","--threads",action="store",help="Number of threads")
    parser.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")
    parser.add_argument("-pv","--pangoLEARN-version", action='version', version=f"pangoLEARN {pangoLEARN.__version__}",help="show pangoLEARN's version number and exit")
    parser.add_argument("--update", action='store_true', default=False, help="Automatically updates to latest release of pangolin and pangoLEARN, then exits")

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)
    args = parser.parse_args()

    if args.update:
        update(__version__, pangoLEARN.__version__)

    snakefile = os.path.join(thisdir, 'scripts','pangolearn.smk')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)
    else:
        print(pfunk.green("Found the snakefile"))

    pfunk.check_installs()

    # to enable not having to pass a query if running update
    # by allowing query to accept 0 to many arguments
    if len(args.query) > 1:
        print(pfunk.cyan(f"Error: Too many query (input) fasta files supplied: {args.query}\nPlease supply one only"))
        parser.print_help()
        sys.exit(-1)
    else:
        # find the query fasta
        query = os.path.join(cwd, args.query[0])
        if not os.path.exists(query):
            sys.stderr.write('Error: cannot find query (input) fasta file at {}\nPlease enter your fasta sequence file and refer to pangolin usage at:\nhttps://github.com/hCoV-2019/pangolin#usage\n for detailed instructions\n'.format(query))
            sys.exit(-1)
        else:
            print(pfunk.green(f"The query file is:") + f"{query}")

        # default output dir
    outdir = ''
    if args.outdir:
        outdir = os.path.join(cwd, args.outdir)
        if not os.path.exists(outdir):
            try:
                os.mkdir(outdir)
            except:
                sys.stderr.write(pfunk.cyan(f'Error: cannot create directory:') + f"{outdir}")
                sys.exit(-1)
    else:
        outdir = cwd


    outfile = ""
    if args.outfile:
        outfile = os.path.join(outdir, args.outfile)
    else:
        outfile = os.path.join(outdir, "lineage_report.csv")

    tempdir = ''
    if args.tempdir:
        to_be_dir = os.path.join(cwd, args.tempdir)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name
    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name

    if args.no_temp:
        print(pfunk.green(f"--no-temp:") + "all intermediate files will be written to {outdir}")
        tempdir = outdir

    if args.alignment:
        align_dir = outdir
        alignment_out = True
    else:
        align_dir = tempdir
        alignment_out = False


    if args.threads:
        print(pfunk.cyan(f"\n--threads flag used, but threading not currently supported. Continuing with one thread."))

    """
    QC steps:
    1) check no empty seqs
    2) check N content
    3) write a file that contains just the seqs to run
    """

    do_not_run = []
    run = []
    for record in SeqIO.parse(query, "fasta"):
        # replace spaces in sequence headers with underscores
        record.id = record.description.replace(' ', '_')
        if "," in record.id:
            record.id=record.id.replace(",","_")

        if len(record) <args.minlen:
            record.description = record.description + f" fail=seq_len:{len(record)}"
            do_not_run.append(record)
            print(record.id, "\tsequence too short")
        else:
            num_N = str(record.seq).upper().count("N")
            prop_N = round((num_N)/len(record.seq), 2)
            if prop_N > args.maxambig:
                record.description = record.description + f" fail=N_content:{prop_N}"
                do_not_run.append(record)
                print(f"{record.id}\thas an N content of {prop_N}")
            else:
                run.append(record)

    if run == []:
        with open(outfile, "w") as fw:
            fw.write("taxon,lineage,probability,pangoLEARN_version,status,note\n")
            for record in do_not_run:
                desc = record.description.split(" ")
                reason = ""
                for item in desc:
                    if item.startswith("fail="):
                        reason = item.split("=")[1]
                fw.write(f"{record.id},None,0,{pangoLEARN.__version__},fail,{reason}\n")
        print(pfunk.cyan(f'Note: no query sequences have passed the qc\n'))
        sys.exit(0)

    post_qc_query = os.path.join(tempdir, 'query.post_qc.fasta')
    with open(post_qc_query,"w") as fw:
        SeqIO.write(run, fw, "fasta")
    qc_fail = os.path.join(tempdir,'query.failed_qc.fasta')
    with open(qc_fail,"w") as fw:
        SeqIO.write(do_not_run, fw, "fasta")

    config = {
        "query_fasta":post_qc_query,
        "outdir":outdir,
        "outfile":outfile,
        "tempdir":tempdir,
        "aligndir":align_dir,
        "alignment_out": alignment_out,
        "trim_start":265,   # where to pad to using datafunk
        "trim_end":29674,   # where to pad after using datafunk
        "qc_fail":qc_fail,
        "pangoLEARN_version":pangoLEARN.__version__
        }

    # find the data
    data_dir = ""
    if args.datadir:
        data_dir = os.path.join(cwd, args.datadir)
        version = "Unknown"
        for r,d,f in os.walk(data_dir):
            for fn in f:
                if fn == "__init__.py":
                    print("Found __init__.py")
                    with open(os.path.join(r, fn),"r") as fr:
                        for l in fr:
                            if l.startswith("__version__"):
                                l = l.rstrip("\n")
                                version = l.split('=')[1]
                                version = version.replace('"',"").replace(" ","")
                                print("pangoLEARN version",version)
        config["pangoLEARN_version"] = version

    if not args.datadir:
        pangoLEARN_dir = pangoLEARN.__path__[0]
        data_dir = os.path.join(pangoLEARN_dir,"data")
    print(f"Looking in {data_dir} for data files...")
    trained_model = ""
    header_file = ""
    lineages_csv = ""

    for r,d,f in os.walk(data_dir):
        for fn in f:
            if fn == "decisionTreeHeaders_v1.joblib":
                header_file = os.path.join(r, fn)
            elif fn == "decisionTree_v1.joblib":
                trained_model = os.path.join(r, fn)
            elif fn == "lineages.metadata.csv":
                lineages_csv = os.path.join(r, fn)
    if trained_model=="" or header_file==""  or lineages_csv=="":
        print(pfunk.cyan("""Check your environment, didn't find appropriate files from the pangoLEARN repo.\n Trained model must be installed, please see https://cov-lineages.org/pangolin.html for installation instructions."""))
        exit(1)
    else:
        if args.decompress:
            prev_size = os.path.getsize(trained_model)

            print("Decompressing model and header files")
            model = joblib.load(trained_model)
            joblib.dump(model, trained_model, compress=0)
            headers = joblib.load(header_file)
            joblib.dump(headers, header_file, compress=0)

            if os.path.getsize(trained_model) >= prev_size:
                print(pfunk.green(f'Success! Decompressed the model file. Exiting\n'))
                sys.exit(0)
            else:
                print(pfunk.cyan(f'Error: failed to decompress model. Exiting\n'))
                sys.exit(0)

        print(pfunk.green("\nData files found"))
        print(f"Trained model:\t{trained_model}")
        print(f"Header file:\t{header_file}")
        print(f"Lineages csv:\t{lineages_csv}")
        config["trained_model"] = trained_model
        config["header_file"] = header_file

    reference_fasta = pkg_resources.resource_filename('pangolin', 'data/reference.fasta')
    config["reference_fasta"] = reference_fasta

    variants_file = pkg_resources.resource_filename('pangolin', 'data/config_b.1.1.7.csv')
    config["b117_variants"] = variants_file

    variants_file = pkg_resources.resource_filename('pangolin', 'data/config_b.1.351.csv')
    config["b1351_variants"] = variants_file

    variants_file = pkg_resources.resource_filename('pangolin', 'data/config_p.1.csv')
    config["p1_variants"] = variants_file

    variants_file = pkg_resources.resource_filename('pangolin', 'data/config_p.2.csv')
    config["p2_variants"] = variants_file

    variants_file = pkg_resources.resource_filename('pangolin', 'data/config_p.3.csv')
    config["p3_variants"] = variants_file


    if args.panGUIlin:
        config["lineages_csv"]=lineages_csv


    if args.verbose:
        quiet_mode = False
        config["log_string"] = ""
    else:
        quiet_mode = True
        lh_path = os.path.realpath(lh.__file__)
        config["log_string"] = f"--quiet --log-handler-script {lh_path} "

    if args.verbose:
        print(pfunk.green("\n**** CONFIG ****"))
        for k in sorted(config):
            print(pfunk.green(k), config[k])

        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                        workdir=tempdir,config=config, cores=1,lock=False
                                        )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=tempdir,
                                    config=config, cores=1,lock=False,quiet=True,log_handler=logger.log_handler
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1


def update(pangolin_version, pangoLEARN_version):
    """
    Using the github releases API check for the latest current release
    of each pangolin and pangoLEARN

    Compare these to the currently running versions and if newer releases
    exist update to them accordingly (or do nothing if current).
    Afterwards, exit program safely with a 0 exit code.

    pangolin_version: string containing the __version__ data for the currently
                      running pangolin module
    pangoLEARN_version: string containing the __version__ data for the imported
                       pangoLEARN data module
    """
    # flag if any element is update if everything is the latest release
    # we want to just continue running
    for dependency, version in [('pangolin', pangolin_version),
                                ('pangoLEARN', pangoLEARN_version)]:
        latest_release = request.urlopen(\
            f"https://api.github.com/repos/cov-lineages/{dependency}/releases")
        latest_release = json.load(latest_release)
        latest_release = LooseVersion(latest_release[0]['tag_name'])

        # to match the tag names add a v to the pangolin internal version
        if dependency == 'pangolin':
            version = "v" + version
        # to match the tag names for pangoLEARN add data release
        elif dependency == 'pangoLEARN':
            version = version.replace(' ', ' data release ')

        # convert to LooseVersion to have proper ordering of versions
        # this prevents someone using the latest commit/HEAD from being
        # downgraded to the last stable release
        version = LooseVersion(version)

        if version < latest_release:
            subprocess.run([sys.executable, '-m', 'pip', 'install', '--upgrade',
                            f"git+https://github.com/cov-lineages/{dependency}.git@{latest_release}"],
                            check=True,
                            stdout=subprocess.DEVNULL,
                            stderr=subprocess.DEVNULL)
            print(f"{dependency} updated to {latest_release}", file=sys.stderr)
        elif version > latest_release:
            print(f"{dependency} ({version}) is newer than latest stable "
                  f"release ({latest_release}), not updating.", file=sys.stderr)
        else:
            print(f"{dependency} already latest release ({latest_release})",
                    file=sys.stderr)

    sys.exit(0)

if __name__ == '__main__':
    main()
