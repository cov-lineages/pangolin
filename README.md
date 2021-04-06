# pangolin

![pangolin test output](https://github.com/pvanheus/pangolin/workflows/pangolin/badge.svg)

**P**hylogenetic **A**ssignment of **N**amed **G**lobal **O**utbreak **LIN**eages

<img src="https://github.com/cov-lineages/pangolin/blob/master/docs/logo.png" width="300">

<strong>Full pangolin documentation found at [cov-lineages.org](https://cov-lineages.org/pangolin.html)</strong>

<strong>Find the pangolin web application [here](https://pangolin.cog-uk.io/), thanks to the Centre for Genomic Pathogen and Surveillance!</strong>

## Quick links
  * [Full documentation](https://cov-lineages.org/pangolin.html)
  * [Requirements](#requirements)
  * [Install pangolin](#install-pangolin)
  * [Check the install worked](#check-the-install-worked)
  * [Updating pangolin](#updating-pangolin)
  * [Updating from pangolin v1.0 to pangolin v2.0](#updating-from-pangolin-v10-to-pangolin-v20)
  * [Basic usage](#basic-usage)
  * [Output](#output)
  * [pangoLEARN description](#pangolearn-description)
  * [Citing pangolin](#citing-pangolin)
  * [References](#references)

### Requirements

Pangolin runs on MacOS and Linux. The conda environment recipe may not build on Windows (I haven't tested it) but can be run using the Windows subsystem for Linux.

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Your query fasta file

### Install pangolin

1. Clone this repository and ``cd pangolin``
2. ``conda env create -f environment.yml``
3. ``conda activate pangolin``
4. ``python setup.py install``
5. That's it

> Troubleshooting install see the [pangolin wiki](https://github.com/cov-lineages/pangolin/wiki)

> Note: we recommend using pangolin in the conda environment specified in the ``environment.yml`` file as per the instructions above. If you can't use conda for some reason, bear in mind the data files are  hosted in two separate repositories at
- [<strong>cov-lineages/lineages</strong>](https://github.com/cov-lineages/lineages.git) 
- [<strong>cov-lineages/pangoLEARN</strong>](https://github.com/cov-lineages/pangoLEARN.git) <br>
you will need to pip install them alongside the other dependencies for pangolin (details found in [<strong>environment.yml</strong>](https://github.com/cov-lineages/pangolin/blob/master/environment.yml)).

### Check the install worked

Type (in the <strong>pangolin</strong> environment):

```
pangolin -v
pangolin -pv
```
and you should see the versions of <strong>pangolin</strong>, and <strong>pangoLEARN</strong> data release printed respectively.


### Updating pangolin

> Note: Even if you have previously installed <strong>pangolin</strong>, as it is being worked on intensively, we recommend you check for updates before running.

To update pangolin and pangoLEARN automatically to the latest stable release:

1. ``conda activate pangolin``
2. ``pangolin --update``

If extra dependencies are introduced (for major releases) the full environment will need to be updated as below:

Alternatively, this can be done manually:

1. ``conda activate pangolin``
2. ``git pull`` \
pulls the latest changes from github
3. ``python setup.py install`` \
re-installs pangolin.
4. ``conda env update -f environment.yml`` \
updates the conda environment (you're unlikely to need to do this, but just in case!)
5. ``pip install git+https://github.com/cov-lineages/pangoLEARN.git --upgrade`` \
updates if there is a new data release

### Updating from pangolin v1.0 to pangolin v2.0 

1. If invoking data path (-d), changed to pangoLEARN instead of lineages
```
-d /home/vix/miniconda3/envs/pangolin/lib/python3.6/site-packages/pangoLEARN/data
```
2. The columns in the output file has also changed, unless running `--legacy`
  - No longer `UFBootstrap`, `aLRT` or `lineages_version`
  - New fields: `probability` and `pangoLEARN_version`

### Basic usage

1. Activate the environment ``conda activate pangolin``
2. Run ``pangolin <query>``, where ``<query>`` is the name of your input file.


### Output

Your output will be a csv file with taxon name and lineage assigned, one line corresponding to each sequence in the fasta file provided

Example:

| Taxon       | Lineage   | support | pangoLEARN_version |status | note |
| ----------- |:---------:|:----------:| :----------:|:----------:| :----------:|
| Virus1      |  B.1      | 80       | 2020-04-27 | passed_qc    | |
| Virus2      |  A.1      |  65      | 2020-04-27 | passed_qc    | |
| Virus3      |  A.3      |  100     | 2020-04-27 | passed_qc    | |
| Virus4      |  B.1.4    |  82      | 2020-04-27 | passed_qc    | |
| Virus5      | None      | 0       | 2020-04-27 | fail    |   N_content:0.80 |
| Virus6      | None      | 0       | 2020-04-27 | fail    |   seq_len:0 |
| Virus7      | None      | 0       | 2020-04-27 | fail    |   failed to map |


### Citing <strong>pangolin</strong>

There is a publication in prep for <strong>pangolin</strong>, but in the meantime please to link to this github [github.com/cov-lineages/pangolin](github.com/cov-lineages/pangolin) if you have used <strong>pangolin</strong> in your research. 

### References

The following external software is run as part of pangolin:

[minimap2](https://github.com/lh3/minimap2)

Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191

[snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.
