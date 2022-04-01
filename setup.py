from setuptools import setup, find_packages
import glob
import os
import pkg_resources

from pangolin import __version__, _program

setup(name='pangolin',
      version=__version__,
      packages=find_packages(),
      scripts=[
            'pangolin/scripts/pangolearn.smk',
            'pangolin/scripts/usher.smk',
            'pangolin/scripts/preprocessing.smk'
                ],
      package_data={"pangolin":["data/*"]},
      install_requires=[
            "biopython>=1.70",
            'pandas>=1.0.1',
            "wheel>=0.34",
            'joblib>=0.11',
            'scikit-learn>=0.23.1',
            "PuLP>=2"
        ],
      description='phylogenetic assignment of named global outbreak lineages',
      url='https://github.com/cov-lineages/pangolin',
      author='Aine OToole & Emily Scher',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = pangolin.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)
