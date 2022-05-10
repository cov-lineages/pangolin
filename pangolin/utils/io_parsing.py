import sys
import os
from pangolin.utils.log_colours import green,cyan
import select
from Bio import SeqIO
import gzip
import lzma

import tempfile
import shutil

from pangolin.utils.config import *


def find_query_file(cwd, tempdir, query_arg):
    if len(query_arg) > 1:
        print(cyan(f"Error: Too many query (input) fasta files supplied: {query_arg}\nPlease supply one only."))
        sys.exit(-1)

    # find the query fasta
    try:
        if not os.path.exists(os.path.join(cwd, query_arg[0])):
            if select.select([sys.stdin,],[],[],0.0)[0]:
                query = os.path.join(tempdir, "stdin_query.fasta")
                with open(query,"w") as fw:
                    for l in sys.stdin:
                        l= l.rstrip("\n")
                        fw.write(l + '\n')
                
                print(green("Query:\t") + "reading from stdin.")
            elif not select.select([sys.stdin,],[],[],0.0)[0]:
                tried_path = os.path.join(cwd, query_arg[0])
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
            query = os.path.join(cwd, query_arg[0])
            print(green(f"Query file:\t") + f"{query}")
    except IndexError:
        sys.stderr.write(cyan(
            f'Error: input query fasta could not be detected from a filepath or through stdin.\n' +
            'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
            ' for detailed instructions.\n'))
        sys.exit(-1)

    return query


def quick_check_query_file(cwd, query_arg, query):
    input_compression_type = "plaintext"
    if os.path.exists(os.path.join(cwd, query_arg[0])):
        file_ending = query.split(".")[-1]
        if file_ending in ["gz","gzip","tgz"]:
            input_compression_type = "gz"
            query = gzip.open(query, 'rt')
        elif file_ending in ["xz","lzma"]:
            input_compression_type = "xz"
            query = lzma.open(query, 'rt')
    try:
        parse= True
        c = 0
        
        for record in SeqIO.parse(query, "fasta"):
            if parse == False:
                break
            parse = False

        return input_compression_type
    except UnicodeDecodeError:
        sys.stderr.write(cyan(
            f'Error: the input query fasta could not be parsed.\n' +
            'Double check your query fasta and that compressed stdin was not passed.\n' +
            'Please enter your fasta sequence file and refer to pangolin usage at: https://cov-lineages.org/pangolin.html' +
            ' for detailed instructions.\n'))
        sys.exit(-1)

def set_up_outdir(outdir_arg,cwd,outdir):
    if outdir_arg:
        outdir = os.path.join(cwd, outdir_arg)
        if not os.path.exists(outdir):
            try:
                os.mkdir(outdir)
            except:
                sys.stderr.write(cyan(f'Error: cannot create directory:') + f"{outdir}")
                sys.exit(-1)
    return outdir

def set_up_outfile(outfile_arg, outfile, outdir):
    if outfile_arg:
        outfile = os.path.join(outdir, outfile_arg)
    else:
        outfile = os.path.join(outdir, outfile)
    return outfile


def set_up_tempdir(tempdir_arg,no_temp_arg,cwd,outdir,config):

    if no_temp_arg:
        tempdir = outdir
        config[KEY_TEMPDIR] = tempdir
        print(green(f"\n--no-temp: ") + f"all intermediate files will be written to {outdir}\n")
    elif tempdir_arg:
        to_be_dir = os.path.join(cwd, tempdir_arg)
        try:
            if not os.path.exists(to_be_dir):
                os.mkdir(to_be_dir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {to_be_dir}.\n'))
            sys.exit(-1)
        tempdir = tempfile.mkdtemp(dir=to_be_dir)
        config[KEY_TEMPDIR] = tempdir
    else:
        tempdir = tempfile.mkdtemp()
        config[KEY_TEMPDIR] = tempdir
        try:
            if not os.path.exists(tempdir):
                os.mkdir(tempdir)
        except:
            sys.stderr.write(cyan(f'Error: cannot create temp directory {tempdir}.\n'))
            sys.exit(-1)
        
    try:
        # write a minimal "constellations" module that scorpio will pick up
        # and from which it will be able to discover constellation files;
        # at the same time, this serves as a test that the tempdir is writable
        constellations_hook = os.path.join(config[KEY_TEMPDIR], 'constellations.py')
        with open(constellations_hook, 'w') as entry_point:
            entry_point.write(f"__version__ = '{config[KEY_CONSTELLATIONS_VERSION]}'\n")
            entry_point.write(f"__path__ = ['{config[KEY_CONSTELLATIONS_PATH]}']\n")
    except:
        sys.stderr.write(cyan(f'Error: cannot write to temp directory {tempdir}.\n'))
        sys.exit(-1)

def cleanup(no_temp,tempdir):
    if not no_temp:
        shutil.rmtree(tempdir)

def parse_alignment_options(alignment_arg, outdir, tempdir,alignment_file_arg, alignment_file):
    if alignment_arg:
        aligndir = outdir
        alignment_out = True
    else:
        aligndir = tempdir
        alignment_out = False

    if alignment_file_arg:
        alignment_file = os.path.join(aligndir, alignment_file_arg)
    else:
        alignment_file = os.path.join(aligndir, alignment_file)

    return alignment_file,alignment_out
