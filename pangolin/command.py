#!/usr/bin/env python3
from pangolin import __version__
import os
import sys
import argparse
import itertools
import snakemake
from tempfile import TemporaryDirectory, TemporaryFile, gettempdir, tempdir
import tempfile
import gzip
import joblib
from pangolin.utils.log_colours import green,cyan,red
import select

try:
    import pangoLEARN
except:
    sys.stderr.write(cyan('Error: please install `pangoLEARN` with \n') +
    "pip install git+https://github.com/cov-lineages/pangoLEARN.git")
    sys.exit(-1)

try:
    import scorpio
except:
    sys.stderr.write(cyan('Error: please install `scorpio` with \n') +
    "pip install git+https://github.com/cov-lineages/scorpio.git")
    sys.exit(-1)

try:
    from pangoLEARN import PANGO_VERSION
except:
    sys.stderr.write(cyan('Error: please update to pangoLEARN version >= 2021-05-27\n'))
    sys.exit(-1)

try:
    import constellations
except:
    sys.stderr.write(cyan('Error: please install `constellations` with \n') +
    "pip install git+https://github.com/cov-lineages/constellations.git")
    sys.exit(-1)

try:
    import pango_designation
except:
    sys.stderr.write(cyan('Error: please install `pango_designation` with \n') +
    "pip install git+https://github.com/cov-lineages/pango-designation.git")
    sys.exit(-1)

from pangolin.utils import dependency_checks
from pangolin.utils import data_install_checks
from pangolin.utils import update

import pangolin.utils.custom_logger as custom_logger

import pkg_resources
from Bio import SeqIO

from . import _program

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def version_from_init(init_file):
    with open(init_file, "r") as fr:
        for l in fr:
            if l.startswith("__version__"):
                l = l.rstrip("\n")
                version = l.split('=')[1]
                version = version.replace('"',"").replace(" ","")
                break
    return version

def main(sysargs = sys.argv[1:]):
    parser = argparse.ArgumentParser(prog = _program,
    description='pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages',
    usage='''pangolin <query> [options]''')

    parser.add_argument('query', nargs="*", help='Query fasta file of sequences to analyse.')
    parser.add_argument('--alignment', action="store_true",help="Optional alignment output.")
    parser.add_argument('--usher', action="store_true",help="Use UShER model instead of default pangoLEARN")
    parser.add_argument('--usher-tree', action='store', dest='usher_protobuf', help="UShER Mutation Annotated Tree protobuf file to use instead of --usher default from pangoLEARN repository or --datadir")
    parser.add_argument('--max-ambig', action="store", default=0.3, type=float,help="Maximum proportion of Ns allowed for pangolin to attempt assignment. Default: 0.3",dest="maxambig")
    parser.add_argument('--min-length', action="store", default=25000, type=int,help="Minimum query length allowed for pangolin to attempt assignment. Default: 25000",dest="minlen")
    parser.add_argument('-o','--outdir', action="store",help="Output directory. Default: current working directory")
    parser.add_argument('--outfile', action="store",help="Optional output file name. Default: lineage_report.csv")
    parser.add_argument('--tempdir',action="store",help="Specify where you want the temp stuff to go. Default: $TMPDIR")
    parser.add_argument("--no-temp",action="store_true",help="Output all intermediate files, for dev purposes.")
    parser.add_argument('-d', '--datadir', action='store',dest="datadir",help="Data directory minimally containing a fasta alignment and guide tree")
    parser.add_argument('--decompress-model',action="store_true",dest="decompress",help="Permanently decompress the model file to save time running pangolin.")
    parser.add_argument("--verbose",action="store_true",help="Print lots of stuff to screen")
    parser.add_argument("-t","--threads",action="store",default=1,type=int, help="Number of threads")
    parser.add_argument("-v","--version", action='version', version=f"pangolin {__version__}")
    parser.add_argument("-pv","--pangoLEARN-version", action='version', version=f"pangoLEARN {pangoLEARN.__version__}",help="show pangoLEARN's version number and exit")
    parser.add_argument("-dv","--pango-designation-version", action='version', version=f"pango-designation {PANGO_VERSION} used for pangoLEARN and UShER training, alias version {pango_designation.__version__}",help="show pango-designation version number used for training and aliases, then exit")
    parser.add_argument("--aliases", action='store_true', default=False, help="print pango-designation alias_key.json and exit")
    parser.add_argument("--skip-designation-hash", action='store_true', default=False, help="Developer option - do not use designation hash to assign lineages")
    parser.add_argument("--update", action='store_true', default=False, help="Automatically updates to latest release of pangolin, pangoLEARN and constellations, then exits")
    parser.add_argument("--update-data", action='store_true',dest="update_data", default=False, help="Automatically updates to latest release of pangoLEARN and constellations, then exits")
    parser.add_argument("--all-versions", action='store_true',dest="all_versions", default=False, help="Print all tool, dependency, and data versions then exit.")

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    # find the data
    if args.datadir is not None:
        # this needs to be an absolute path when we pass it to scorpio
        args.datadir = os.path.abspath(args.datadir)
    alias_file = None
    pango_designation_dir = pango_designation.__path__[0]
    constellations_dir = constellations.__path__[0]
    constellation_files = []
    data_locations = [os.walk(pango_designation_dir), os.walk(constellations_dir)]
    if args.datadir is not None:
        data_locations.append(os.walk(args.datadir))
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
                pango_designation.__version__ = version_from_init(os.path.join(r, '__init__.py'))
            elif r.endswith('/constellations') and fn == '__init__.py':
                constellations.__version__ = version_from_init(os.path.join(r, fn))
            elif (r.endswith('/constellations') or r.endswith('/constellations/definitions')) and fn.endswith('.json'):
                constellation_files.append(os.path.join(r, fn))


    use_datadir = False
    if args.datadir:
        data_dir = os.path.join(cwd, args.datadir)
        version = "Unknown"
        for r,d,f in os.walk(data_dir):
            for fn in f:
                if r.endswith('pangoLEARN') and fn == "__init__.py":
                    # print("Found __init__.py")
                    version = version_from_init(os.path.join(r, fn))
                    if version > pangoLEARN.__version__:
                        # only use this for pangoLEARN if the version is > than what we already have
                        pangoLEARN.__version__ = version
                        use_datadir = True
    if use_datadir == False:
        # we haven't got a viable datadir from searching args.datadir
        pangoLEARN_dir = pangoLEARN.__path__[0]
        data_dir = os.path.join(pangoLEARN_dir,"data")

    # print(f"Looking in {data_dir} for data files...")
    if args.update:
        update.update({'pangolin': __version__,
                'pangolearn': pangoLEARN.__version__,
                'constellations': constellations.__version__,
                'scorpio': scorpio.__version__,
                'pango-designation': pango_designation.__version__
                })

    if args.update_data:
        update.update({'pangolearn': pangoLEARN.__version__,
                'constellations': constellations.__version__,
                'pango-designation': pango_designation.__version__}, args.datadir)

    if not alias_file:
        sys.stderr.write(cyan('Could not find alias file: please update pango-designation with \n') +
                         "pip install git+https://github.com/cov-lineages/pango-designation.git")
        sys.exit(-1)

    if args.aliases:
        with open(alias_file, 'r') as handle:
            for line in handle:
                print(line.rstrip())
        sys.exit(0)

    if args.all_versions:
        print(f"pangolin: {__version__}\n"
              f"pangolearn: {pangoLEARN.__version__}\n"
              f"constellations: {constellations.__version__}\n"
              f"scorpio: {scorpio.__version__}\n"
              f"pango-designation used by pangoLEARN/Usher: {PANGO_VERSION}\n"
              f"pango-designation aliases: {pango_designation.__version__}")
        sys.exit(0)


    dependency_checks.check_dependencies(args.usher)

    # to enable not having to pass a query if running update
    # by allowing query to accept 0 to many arguments
    if len(args.query) > 1:
        print(cyan(f"Error: Too many query (input) fasta files supplied: {args.query}\nPlease supply one only"))
        parser.print_help()
        sys.exit(-1)
    else:
        # find the query fasta
        if not args.decompress:
            try:
                if not os.path.exists(os.path.join(cwd, args.query[0])):
                    if select.select([sys.stdin,],[],[],0.0)[0]:
                        query = sys.stdin
                    elif not select.select([sys.stdin,],[],[],0.0)[0]:
                        tried_path = os.path.join(cwd, args.query[0])
                        if tried_path.endswith("-"):
                            sys.stderr.write(cyan(
                                f'Error: cannot find query (input) fasta file using stdin.\n' +
                                         'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
                                         ' for detailed instructions.\n'))
                            sys.exit(-1)
                        else:
                            sys.stderr.write(cyan(f'Error: cannot find query (input) fasta file at:') + f'{tried_path}\n' +
                                          'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
                                          ' for detailed instructions.\n')
                            sys.exit(-1)
                else:
                    query = os.path.join(cwd, args.query[0])
                    print(green(f"The query file is:") + f"{query}")
            except IndexError:
                sys.stderr.write(cyan(
                    f'Error: input query fasta could not be detected from a filepath or through stdin.\n' +
                    'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
                    ' for detailed instructions.\n'))
                sys.exit(-1)

        # default output dir

    if args.outdir:
        outdir = os.path.join(cwd, args.outdir)
        if not os.path.exists(outdir):
            try:
                os.mkdir(outdir)
            except:
                sys.stderr.write(cyan(f'Error: cannot create directory:') + f"{outdir}")
                sys.exit(-1)
    else:
        outdir = cwd

    if args.outfile:
        outfile = os.path.join(outdir, args.outfile)
    else:
        outfile = os.path.join(outdir, "lineage_report.csv")

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
        print(green(f"\n--no-temp: ") + f"all intermediate files will be written to {outdir}\n")
        tempdir = outdir

    if args.alignment:
        align_dir = outdir
        alignment_out = True
    else:
        align_dir = tempdir
        alignment_out = False

    """
    QC steps:
    1) check no empty seqs
    2) check N content
    3) write a file that contains just the seqs to run
    """
    if not args.decompress:
        do_not_run = []
        run = []
        
        print(green("** Running sequence QC **"))

        if not select.select([sys.stdin,],[],[],0.0)[0]:
            file_ending = query.split(".")[-1]
            if file_ending in ["gz","gzip","tgz"]:
                query = gzip.open(query, 'rt')
            elif file_ending in ["xz","lzma"]:
                query = lzma.open(query, 'rt')
                
        post_qc_query = os.path.join(tempdir, 'query.post_qc.fasta')
        fw_pass = open(post_qc_query,"w")
        qc_fail = os.path.join(tempdir,'query.failed_qc.fasta')
        fw_fail = open(qc_fail,"w")

        total_input = 0
        total_pass = 0
        
        try:
            for record in SeqIO.parse(query, "fasta"):
                total_input +=1
                record.description = record.description.replace(' ', '_').replace(",","_")
                record.id = record.description
                if "," in record.id:
                    record.id=record.id.replace(",","_")

                if len(record) <args.minlen:
                    record.description = record.description + f" fail=seq_len:{len(record)}"
                    fw_fail.write(f">{record.description}\n{record.seq}\n")
                else:
                    num_N = str(record.seq).upper().count("N")
                    prop_N = round((num_N)/len(record.seq), 2)
                    if prop_N > args.maxambig:
                        record.description = record.description + f" fail=N_content:{prop_N}"
                        fw_fail.write(f">{record.description}\n{record.seq}\n")
                    else:
                        total_pass +=1
                        seq = str(record.seq).replace("-","")
                        fw_pass.write(f">{record.description}\n{seq}\n")
        except UnicodeDecodeError:
            sys.stderr.write(cyan(
                f'Error: input query fasta could not be detected from a filepath or through stdin.\n' +
                'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
                ' for detailed instructions.\n'))
            sys.exit(-1)

        print(green("Number of sequences detected: ") + f"{total_input}")
        print(green("Total passing QC: ") + f"{total_pass}")
        fw_fail.close()
        fw_pass.close()

        if total_pass == 0:
            with open(outfile, "w") as fw:
                fw.write("taxon,lineage,conflict,ambiguity_score,scorpio_call,scorpio_support,scorpio_conflict,version,pangolin_version,pangoLEARN_version,pango_version,status,note\n")
                for record in SeqIO.parse(os.path.join(tempdir,'query.failed_qc.fasta'), "fasta"):
                    desc = record.description.split(" ")
                    reason = ""
                    for item in desc:
                        if item.startswith("fail="):
                            reason = item.split("=")[1]
                    fw.write(f"{record.id},None,,,,,,PANGO-{PANGO_VERSION},{__version__},{pangoLEARN.__version__},{PANGO_VERSION},fail,{reason}\n")
            print(cyan(f'Note: no query sequences have passed the qc\n'))
            sys.exit(0)

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
            "alias_file": alias_file,
            "constellation_files": constellation_files,
            "skip_designation_hash": args.skip_designation_hash,
            "verbose":args.verbose,
            "pangoLEARN_version":pangoLEARN.__version__,
            "pangolin_version":__version__,
            "pango_version":PANGO_VERSION,
            "threads":args.threads
            }

        data_install_checks.check_install(config)
        snakefile = data_install_checks.get_snakefile(thisdir)

        dependency_checks.set_up_verbosity(config)

    trained_model = ""
    header_file = ""
    designated_hash=""
    use_usher = args.usher
    if args.usher_protobuf:
        usher_protobuf = os.path.join(cwd, args.usher_protobuf)
        if not os.path.exists(usher_protobuf):
            sys.stderr.write('Error: cannot find --usher-tree file at {}\n'.format(usher_protobuf))
            sys.exit(-1)
        use_usher = True
    else:
        usher_protobuf = ""

    for r,d,f in os.walk(data_dir):
        for fn in f:
            if fn == "decisionTreeHeaders_v1.joblib":
                header_file = os.path.join(r, fn)
            elif fn == "decisionTree_v1.joblib":
                trained_model = os.path.join(r, fn)
            elif fn =="lineages.hash.csv":
                designated_hash = os.path.join(r, fn)
            elif fn == "lineageTree.pb" and usher_protobuf == "":
                usher_protobuf = os.path.join(r, fn)
    if ((use_usher and (usher_protobuf == "" or designated_hash=="") or
        (not use_usher and (trained_model=="" or header_file=="" or designated_hash=="")))):
        print(cyan("""pangoLEARN version should be >= 2021-05-27. \n
Appropriate data files not found from the installed pangoLEARN repo.
Please see https://cov-lineages.org/pangolin.html for installation and updating instructions."""))
        exit(1)
    else:
        if args.decompress:
            prev_size = os.path.getsize(trained_model)

            print("Decompressing model and header files.")
            model = joblib.load(trained_model)
            joblib.dump(model, trained_model, compress=0)
            headers = joblib.load(header_file)
            joblib.dump(headers, header_file, compress=0)

            if os.path.getsize(trained_model) >= prev_size:
                print(green(f'Success! Decompressed the model file. Exiting\n'))
                sys.exit(0)
            else:
                print(cyan(f'Error: failed to decompress model. Exiting\n'))
                sys.exit(0)

        print(green("\nData files found:"))
        if use_usher:
            print(f"UShER tree:\t{usher_protobuf}")
            print(f"Designated hash:\t{designated_hash}")
        else:
            print(f"Trained model:\t{trained_model}")
            print(f"Header file:\t{header_file}")
            print(f"Designated hash:\t{designated_hash}")

        config["trained_model"] = trained_model
        config["header_file"] = header_file
        config["designated_hash"] = designated_hash

    if use_usher:
        config["usher_protobuf"] = usher_protobuf

    if config['verbose']:
        print(green("\n**** CONFIG ****"))
        for k in sorted(config):
            print(green(k), config[k])

        status = snakemake.snakemake(snakefile, printshellcmds=True, forceall=True, force_incomplete=True,
                                        workdir=tempdir,config=config, cores=args.threads,lock=False
                                        )
    else:
        logger = custom_logger.Logger()
        status = snakemake.snakemake(snakefile, printshellcmds=False, forceall=True,force_incomplete=True,workdir=tempdir,
                                    config=config, cores=args.threads,lock=False,quiet=True,log_handler=config["log_api"]
                                    )

    if status: # translate "success" into shell exit code of 0
       return 0

    return 1


if __name__ == '__main__':
    main()
