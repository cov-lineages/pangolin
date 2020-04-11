from setuptools import setup, find_packages
import glob
import os
import pkg_resources
# Note: the _program variable is set in __init__.py.
# it determines the name of the package/final command line tool.
from pangolin import __version__, _program

setup(name='pangolin',
      version=__version__,
      packages=find_packages(),
      scripts=['pangolin/scripts/assign_query_file.smk',
                'pangolin/scripts/assign_query_lineage.smk',
                'pangolin/scripts/Snakefile',
                'pangolin/scripts/assign_lineage.py',
                'pangolin/scripts/lineage_finder.py',
                'pangolin/scripts/utils.py'
                ],
      package_data={'pangolin':['config.yaml',
                                'data/*']},
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0"
        ],
      description='hcov-2019 subtyping command line tool',
      url='https://github.com/aineniamh/lineages',
      author='Aine OToole and JT McCrone',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = pangolin.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)