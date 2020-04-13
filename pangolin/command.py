#!/usr/bin/env python3

import argparse
import os.path
import snakemake
import sys
import pprint
import json

from . import _program


thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='pangolin: Pipeline for Assigning Global Outbreak Lineages', 
    usage='''pangolin <query> [options]''')

    parser.add_argument('query')
    parser.add_argument('-o','--outdir', action="store")
    parser.add_argument('-n', '--dry-run', action='store_true')
    parser.add_argument('-f', '--force', action='store_true')
    parser.add_argument('-t', '--threads', action='store',type=int)
    parser.add_argument('-u','--unlock', action='store_true')
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

    # run subtyping
    status = snakemake.snakemake(snakefile, printshellcmds=True,
                                 dryrun=args.dry_run, forceall=args.force,
                                 config=config, unlock=args.unlock, cores=threads,lock=False
                                 )

    if status: # translate "success" into shell exit code of 0
       return 0
    return 1

if __name__ == '__main__':
    main()