# pangolin

Phylogenetic Assignment of Named Global Outbreak LINeages

### Requirements

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Your query fasta file

### Install pangolin

1. Clone this repository
2. ``conda env create -f environment.yml``
3. That's it

### Usage

1. Activate the environment ``conda activate pangolin-env``
2. Run ``lineage <query>``

```

pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages

usage: pangolin <query> [options]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR, --outdir OUTDIR
  -n, --dry-run
  -f, --force
  -t THREADS, --threads THREADS
  -u, --unlock
  ```

### Output

Your output will be a csv file with taxon name and lineage assigned, one line corresponding to each sequence in the fasta file provided

Example:

| Taxon       | Lineage   |
| ----------- |:---------:|
| Virus1      |  B.1      |
| Virus2      |  A.1      |
| Virus3      |  A.3      |
| Virus4      |  B.1.4    |
