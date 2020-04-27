#!/usr/bin/env python3
from pangolin import __version__
import argparse
import os.path
import snakemake
import sys
import pprint
import json
from Bio import SeqIO

from . import _program


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages', 
    usage='''pangolin <query> [options]''')

    parser.add_argument('query')
    parser.add_argument('-o','--outdir', action="store",help="Output directory")
    parser.add_argument('-d', '--data', action='store',help="Data directory minimally containing a fasta alignment and guide tree")
    parser.add_argument('-n', '--dry-run', action='store_true',help="Go through the motions but don't actually run")
    parser.add_argument('-f', '--force', action='store_true',help="Overwrite all output")
    parser.add_argument('-t', '--threads', action='store',type=int,help="Number of threads")
    parser.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # find the Snakefile
    snakefile = os.path.join(thisdir, 'scripts/Snakefile')
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

    for record in SeqIO.parse(query, "fasta"):
        num_N = str(record.seq).upper().count("N")
        pcent_N = (num_N*100)/len(record.seq)
        if pcent_N > 50: 
            print(f"Error: {record.id} has an N content greater than 50%, which may lead to inaccurate assignment.\nPlease remove this sequence and try again.\nExiting.",file=sys.stderr)
            exit(1)

    # default output dir
    if args.outdir:
        outdir = args.outdir.rstrip("/")
    else:
        outdir = cwd.rstrip("/")

    # how many threads to pass
    if args.threads:
        threads = args.threads
    else:
        threads = 1

    print("Number of threads is", threads)

    config = {
        "query_fasta":query,
        "outdir":outdir
        }

    if args.data:
        data_dir = os.path.join(cwd, args.data.rstrip("/"))
        print(f"Looking in {data_dir} for data files...")
        representative_aln = ""
        guide_tree = ""
        lineages_csv = ""
        for r,d,f in os.walk(data_dir):
            for fn in f:
                if fn.endswith(".fasta"):
                    representative_aln = r + '/' + fn
                elif fn.endswith(".tree") or fn.endswith(".treefile"):
                    guide_tree = r + '/' + fn
                elif fn.endswith(".csv"):
                    lineages_csv = r + "/" + fn
        print("\nData files found")
        print(f"Sequence alignment:\t{representative_aln}")
        print(f"Guide tree:\t\t{guide_tree}")
        print(f"Lineages csv:\t\t{lineages_csv}")
        if representative_aln=="" or guide_tree=="":
            print("Didn't find appropriate files.\nTreefile must end with `.tree` or `.treefile`.\nAlignment must be in `.fasta` format.\nExiting.")
            exit(1)
        else:
            config["representative_aln"]=representative_aln
            config["guide_tree"]=guide_tree
            config["lineages_csv"]=lineages_csv

    # run subtyping
    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force,force_incomplete=True,
                                 config=config, cores=threads,lock=False,quiet=True
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0
    return 1

if __name__ == '__main__':
    main()