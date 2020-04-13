# pangolin

Phylogenetic Assignment of Named Global Outbreak LINeages

### Requirements

Pangolin runs on MacOS and Linux. The conda environment recipe may not build on Windows (I haven't tested it) but can be run using the Windows subsystem for Linux.

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Your query fasta file

### Install pangolin

1. Clone this repository
2. ``cd pangolin``
2. ``conda env create -f environment.yml``
3. ``python setup.py install``
3. That's it

### Usage

1. Activate the environment ``conda activate pangolin``
2. Run ``pangolin <query>``

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
