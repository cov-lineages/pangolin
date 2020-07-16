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
                'pangolin/scripts/prepare_package_data.smk',
                'pangolin/scripts/Snakefile',
                'pangolin/scripts/assign_lineage.py',
                'pangolin/scripts/lineage_finder.py',
				# 'pangolin/scripts/codify_all_sites.py',
                'pangolin/scripts/utils.py',
				'pangolin/scripts/get_polytomy.py',
				'pangolin/scripts/get_basal_polytomy.py',
                'pangolin/scripts/categorise_snps.py',
                'pangolin/scripts/find_all_snps.py',
                'pangolin/scripts/prepare_package_data.smk',
                'pangolin/scripts/get_masked_representatives.py',
                'pangolin/scripts/sam_2_snps.py',
                'pangolin/scripts/snp_based_classifier.py',
                'pangolin/scripts/snp_based_classify.smk',
                'pangolin/scripts/report_classes.py',
                'pangolin/scripts/report_results.py',
                "pangolin/scripts/get_basal_polytomy.py",
                'pangolin/scripts/get_polytomy.py'
                ],
      install_requires=[
            "biopython>=1.70",
            "dendropy>=4.4.0",
            "pytools>=2020.1",
            'pandas>=1.0.1',
            'pysam>=0.15.4'
        ],
      description='hcov-2019 subtyping command line tool',
      url='https://github.com/hCov-2019/pangolin',
      author='Aine OToole and JT McCrone',
      author_email='aine.otoole@ed.ac.uk',
      entry_points="""
      [console_scripts]
      {program} = pangolin.command:main
      """.format(program = _program),
      include_package_data=True,
      keywords=[],
      zip_safe=False)