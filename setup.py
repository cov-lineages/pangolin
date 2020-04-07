from setuptools import setup, find_packages
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from lineage import __version__, _program

setup(name='lineage',
      version=__version__,
      packages=['lineage'],
      scripts=['lineage/bin/assign_query_file.smk',
                'lineage/bin/assign_query_lineage.smk',
                'lineage/bin/Snakefile'],
      package_data={'lineage':['config.yaml',
                                'data/*']},
      description='hcov-2019 subtyping command line tool',
      url='',
      author='Aine O Toole',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = lineage.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)