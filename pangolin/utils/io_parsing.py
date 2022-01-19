import sys
import os
from pangolin.utils.log_colours import green,cyan
import select


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