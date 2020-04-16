# pangolin

**P**hylogenetic **A**ssignment of **N**amed **G**lobal **O**utbreak **LIN**eages

<img src="https://github.com/hCoV-2019/pangolin/blob/master/docs/logo.png" width="300">

### Requirements

Pangolin runs on MacOS and Linux. The conda environment recipe may not build on Windows (I haven't tested it) but can be run using the Windows subsystem for Linux.

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Your query fasta file

### Install pangolin

1. Clone this repository and ``cd pangolin``
2. ``conda env create -f environment.yml``
3. ``python setup.py install`` or ``pip install .``
4. That's it

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

| Taxon       | Lineage   | UFbootstrap |
| ----------- |:---------:|:----------:|
| Virus1      |  B.1      |  82     |
| Virus2      |  A.1      |  95     |
| Virus3      |  A.3      |  100    |
| Virus4      |  B.1.4    |  73     |


### Authors

Pangolin was created by Áine O'Toole and JT McCrone

