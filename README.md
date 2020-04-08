# lineages
Resources for a lineage naming scheme for SARS-CoV-2/hCoV-2019


## lineage assignment tool

### Requirements

Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html).

### Install lineage assignment tool

1. Clone this repository.
2. ``conda env create -f environment.yml``
3. That's it. *

    \* caveat it's actually not it right now...

4. Download the 2.0-rc2 version of iqtree [here](http://www.iqtree.org/#download).
5. Make sure iqtree 2 lives in your path.

>Recommendation: Replace the downloaded iqtree executable with the one in the lineage-env environment directory. 


### Run lineage assignment tool

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

### Output 

Your output will be a csv file with taxon name and lineage assigned, one line corresponding to each sequence in the fasta file provided. 

Example:

| Taxon        | Lineage      |
| ------------- |:-------------:|
| Virus1      | B.1      |
| Virus2      | A.1      |
| Virus3      | A.3      |
| Virus4      | B.1.4      |
