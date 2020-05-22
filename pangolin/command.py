#!/usr/bin/env python3
from pangolin import __version__
import argparse
import os.path
import snakemake
import sys
from tempfile import gettempdir
import tempfile
import pprint
import json
import lineages
import setuptools
from Bio import SeqIO

from . import _program


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages', 
    usage='''pangolin <query> [options]''')

    parser.add_argument('query')
    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('--outfile', action="store",help="Optional output file name. Default: lineage_report.csv")
    parser.add_argument('-d', '--data', action='store',help="Data directory minimally containing a fasta alignment and guide tree")
    parser.add_argument('-n', '--dry-run', action='store_true',help="Go through the motions but don't actually run")
    parser.add_argument('-f', '--force', action='store_true',help="Overwrite all output",dest="force")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument('--max-ambig', action="store", default=0.5, type=float,help="Maximum proportion of Ns allowed for pangolin to attempt assignment. Default: 0.5",dest="maxambig")
    parser.add_argument('--min-length', action="store", default=10000, type=int,help="Minimum query length allowed for pangolin to attempt assignment. Default: 10000",dest="minlen")
    parser.add_argument('--panGUIlin', action='store_true',help="Run web-app version of pangolin")
    parser.add_argument('--write-tree', action='store_true',help="Output a phylogeny for each query sequence placed in the guide tree",dest="write_tree")
    parser.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    parser.add_argument("-p","--include-putative",action="store_true",help="Include the bleeding edge lineage definitions in assignment",dest="include_putative")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")
    parser.add_argument("-lv","--lineages-version", action='version', version=f"lineages {lineages.__version__}",help="show lineages's version number and exit")

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # find the Snakefile
    snakefile = os.path.join(thisdir, 'scripts','Snakefile')
    if not os.path.exists(snakefile):
        sys.stderr.write('Error: cannot find Snakefile at {}\n'.format(snakefile))
        sys.exit(-1)
    else:
        print("Found the snakefile")

    # find the query fasta
    query = os.path.join(cwd, args.query)
    if not os.path.exists(query):
        sys.stderr.write('Error: cannot find query at {}\n'.format(query))
        sys.exit(-1)
    else:
        print(f"The query file is {query}")

        # default output dir
    outdir = ''
    if args.outdir:
        outdir = os.path.join(cwd, args.outdir)
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

    """ 
    QC steps:
    1) check no empty seqs
    2) check N content
    3) write a file that contains just the seqs to run
    """

    do_not_run = []
    run = []
    for record in SeqIO.parse(query, "fasta"):
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

    post_qc_query = os.path.join(tempdir, 'query.post_qc.fasta')
    with open(post_qc_query,"w") as fw:
        SeqIO.write(run, fw, "fasta")
    qc_fail = os.path.join(tempdir,'query.failed_qc.fasta')
    with open(qc_fail,"w") as fw:
        SeqIO.write(do_not_run, fw, "fasta")

    # how many threads to pass
    if args.threads:
        threads = args.threads
    else:
        threads = 1

    print("Number of threads is", threads)

    config = {
        "query_fasta":post_qc_query,
        "outdir":outdir,
        "outfile":outfile,
        "tempdir":tempdir,
        "qc_fail":qc_fail,
        "lineages_version":lineages.__version__
        }

    if args.force:
        config["force"]="forceall"
    # find the data
    data_dir = ""
    if args.data:
        data_dir = os.path.join(cwd, args.data)
    else:
        lineages_dir = lineages.__path__[0]
        data_dir = os.path.join(lineages_dir,"data")

    print(f"Looking in {data_dir} for data files...")
    representative_aln = ""
    guide_tree = ""
    lineages_csv = ""
    for r,d,f in os.walk(data_dir):
        for fn in f:
            if args.include_putative:
                if fn.endswith("putative.fasta"):
                    representative_aln = os.path.join(r, fn)
                elif fn.endswith("putative.fasta.treefile"):
                    guide_tree = os.path.join(r, fn)
                elif fn.endswith(".csv") and fn.startswith("lineages"):
                    lineages_csv = os.path.join(r, fn)
            else:
                if fn.endswith("safe.fasta"):
                    representative_aln = os.path.join(r, fn)
                elif fn.endswith("safe.fasta.treefile"):
                    guide_tree = os.path.join(r, fn)
                elif fn.endswith(".csv") and fn.startswith("lineages"):
                    lineages_csv = os.path.join(r, fn)

    print("\nData files found")
    print(f"Sequence alignment:\t{representative_aln}")
    print(f"Guide tree:\t\t{guide_tree}")
    print(f"Lineages csv:\t\t{lineages_csv}")
    if representative_aln=="" or guide_tree=="" or lineages_csv=="":
        print("""Didn't find appropriate files.\nTreefile must end with `.treefile`.\nAlignment must be in `.fasta` format.\n \
If you've specified --include-putative \n 
you must have files ending in putative.fasta.treefile\nExiting.""")
        exit(1)
    else:
        config["representative_aln"]=representative_aln
        config["guide_tree"]=guide_tree

    if args.write_tree:
        config["write_tree"]="True"

    if args.panGUIlin:
        config["lineages_csv"]=lineages_csv

    if args.verbose:
        quiet_mode = False
    else:
        quiet_mode = True

    # run subtyping
    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force,force_incomplete=True,
                                 config=config, cores=threads,lock=False,quiet=quiet_mode,workdir=tempdir
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1

if __name__ == '__main__':
    main()