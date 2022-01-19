import sys
from pangolin.utils.log_colours import green,cyan


def install_error(package, url):
    sys.stderr.write(cyan(f'Error: please install `{package}` with \n') +
    f"pip install git+{url}")
    sys.exit(-1)
