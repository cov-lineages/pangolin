import sys
import os
from pangolin.utils.log_colours import green,cyan
import select


from tempfile import TemporaryDirectory, TemporaryFile, gettempdir, tempdir
import tempfile


def find_query_file(cwd, query_arg):
    if len(query_arg) > 1:
        print(cyan(f"Error: Too many query (input) fasta files supplied: {query_arg}\nPlease supply one only."))
        sys.exit(-1)

    # find the query fasta
    try:
        if not os.path.exists(os.path.join(cwd, query_arg[0])):
            if select.select([sys.stdin,],[],[],0.0)[0]:
                query = sys.stdin
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
            print(green(f"The query file is:") + f"{query}")
    except IndexError:
        sys.stderr.write(cyan(
            f'Error: input query fasta could not be detected from a filepath or through stdin.\n' +
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


def set_up_tempdir(tempdir_arg,no_temp_arg,cwd,outdir):

    if tempdir_arg:
        to_be_dir = os.path.join(cwd, tempdir_arg)
        if not os.path.exists(to_be_dir):
            os.mkdir(to_be_dir)
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=to_be_dir)
        tempdir = temporary_directory.name
    else:
        temporary_directory = tempfile.TemporaryDirectory(suffix=None, prefix=None, dir=None)
        tempdir = temporary_directory.name

    if no_temp_arg:
        print(green(f"\n--no-temp: ") + f"all intermediate files will be written to {outdir}\n")
        tempdir = outdir

    return tempdir


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
