# lineages
Resources for a lineage naming scheme for SARS-CoV-2/hCoV-2019

## Install subtyping tool

1. Clone this repository.
2. ``conda env create -f environment.yml``
3. That's it. 

## Run subtyping tool

1. Activate the environment ``conda activate lineage-test``
2. Run ``lineage <query>``
Command line options:
```
lineage: run hcov-2019 lineage assignment tool

positional arguments:
  query

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
  -n, --dry-run
  -f, --force
  -t THREADS, --threads THREADS
  ```