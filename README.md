# pangolin

![pangolin test output](https://github.com/pvanheus/pangolin/workflows/pangolin/badge.svg)
[![BioConda version](https://anaconda.org/bioconda/pangolin/badges/version.svg)](https://anaconda.org/bioconda/pangolin)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=pangolin)


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
