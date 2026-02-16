from setuptools import setup, find_packages

setup(name='pangolin',
      packages=find_packages(),
      scripts=[
            'pangolin/scripts/pangolearn.smk',
            'pangolin/scripts/usher.smk',
            'pangolin/scripts/preprocessing.smk'
      ],
      package_data={"pangolin":["data/*"]}
      )
