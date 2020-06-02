# pangolin

**P**hylogenetic **A**ssignment of **N**amed **G**lobal **O**utbreak **LIN**eages

<img src="https://github.com/hCoV-2019/pangolin/blob/master/docs/logo.png" width="300">


## Quick links

  * [Requirements](#requirements)
  * [Install pangolin](#install-pangolin)
  * [Check the install worked](#check-the-install-worked)
  * [Updating pangolin](#updating-pangolin)
  * [Usage](#usage)
  * [Output](#output)
  * [Recall rate](#recall-rate)
  * [Source data](#source-data)
  * [Authors](#authors)
  * [Citing ``pangolin``](#citing-pangolin)
  * [References](#references)


### Requirements

Pangolin runs on MacOS and Linux. The conda environment recipe may not build on Windows (I haven't tested it) but can be run using the Windows subsystem for Linux.

1. Some version of conda, we use Miniconda3. Can be downloaded from [here](https://docs.conda.io/en/latest/miniconda.html)
2. Your query fasta file

### Install pangolin

1. Clone this repository and ``cd pangolin``
2. ``conda env create -f environment.yml``
3. ``conda activate pangolin``
4. ``python setup.py install`` or ``pip install .``
5. That's it

> Note: we recommend using pangolin in the conda environment specified in the ``environment.yml`` file as per the instructions above. If you can't use conda for some reason, bear in mind the data files are now hosted in a separate repository at [``hCoV-2019/lineages``](https://github.com/hCoV-2019/lineages.git) and you will need to pip install that alongside the other dependencies for pangolin (details found in [``environment.yml``](https://github.com/hCoV-2019/pangolin/blob/master/environment.yml)). 

### Check the install worked

Type (in the pangolin environment):

```
pangolin -v
pangolin -lv
```
and you should see the versions of pangolin and lineages data release printed respectively.

### Updating pangolin

> Note: Even if you have previously installed ``pangolin``, as it is being worked on intensively, we recommend you check for updates before running.

To update:

1. ``conda activate pangolin``
2. ``git pull`` \
pulls the latest changes from github
3. ``python setup.py install`` \
re-installs pangolin
4. ``pip install git+https://github.com/hCoV-2019/lineages.git --upgrade`` \
updates if there is a new data release
5. ``conda env update -f environment.yml`` \
updates the conda environment (you're unlikely to need to do this, but just in case!)


### Usage

1. Activate the environment ``conda activate pangolin``
2. Run ``pangolin <query>``, where ``<query>`` is the name of your input file.

```
pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages

positional arguments:
  input_file - input file with data to analyse (FASTA format)

optional arguments:
  -h, --help                  show this help message and exit
  -o OUTDIR, --outdir OUTDIR  Output directory
  -d DATA, --data DATA        Data directory minimally containing a fasta alignment
                              and guide tree
  -n, --dry-run               Go through the motions but don't actually run
  -f, --force                 Overwrite all output
  --tempdir TEMPDIR           Specify where you want the temp stuff to go. Default:
                              $TMPDIR
  --panGUIlin                 Run web-app version of pangolin
  --max-ambig MAXAMBIG        Maximum proportion of Ns allowed for pangolin to
                              attempt assignment. Default: 0.5
  --min-length MINLEN         Minimum query length allowed for pangolin to attempt
                              assignment. Default: 10000
  -t THREADS, --threads       THREADS
                              Number of threads
  -v, --version               show program's version number and exit
  -lv, --lineages-version     show lineages's version number and exit
  ```

### Output

Your output will be a csv file with taxon name and lineage assigned, one line corresponding to each sequence in the fasta file provided

Example:

| Taxon       | Lineage   | aLRT | UFbootstrap | lineages_version |status | note |
| ----------- |:---------:|:----------:|:----------:| :----------:|:----------:| :----------:|
| Virus1      |  B.1      | 80      |  82    | 2020-04-27 | passed_qc    | |
| Virus2      |  A.1      |  65     | 95     | 2020-04-27 | passed_qc    | |
| Virus3      |  A.3      |  100     | 100    | 2020-04-27 | passed_qc    | |
| Virus4      |  B.1.4    |  82     | 73     | 2020-04-27 | passed_qc    | |
| Virus5      | None      | 0       | 0       | 2020-04-27 | fail    |   N_content:0.80 |
| Virus6      | None      | 0       | 0       | 2020-04-27 | fail    |   seq_len:0 |

Resources for interpreting the aLRT and UFbootstrap output can be found [here](http://www.iqtree.org/doc/Tutorial#assessing-branch-supports-with-single-branch-tests) and [here](http://www.iqtree.org/doc/Command-Reference).

### Recall rate
Of 9,843 GISAID sequences assigned lineages by hand (taking sequence, phylogeny and metadata into account), pangolin accurately assigns the lineage of 97.85% of those sequences. Of the sequences that were not recalled correctly, 74.5% had 0 bootstrap and 0 alrt. We're continuing to work to improve this recall rate, but recommend interpreting the pangolin output cautiously with due attention to the UFbootstrap and aLRT values. 

Given hCoV-2019 is relatively slow evolving for an RNA virus and there is still not a huge amount of diversity, missing or ambiguous data at key residues may lead to incorrect placement within the guide tree. We have a filter in place that by default with not call a lineage for any sequence with >50% N-content, but this can be made more conservative with the command line option `--max-ambig`.

### Source data

``pangolin`` runs using a guide tree and alignment hosted at [``hCoV-2019/lineages``](https://github.com/hCoV-2019/lineages.git). Some of this data is sourced from GISAID, but anonymised and encrypted to fit with guidelines. Appropriate permissions have been given and acknowledgements for the teams that have worked to provide the original SARS-CoV-2 genome sequences to GISAID are also hosted in [``hCoV-2019/lineages``](https://raw.githubusercontent.com/hCoV-2019/lineages/master/gisaid_acknowledgements.tsv).


### Authors

Pangolin was created by [Áine O'Toole](https://aineotoole.co.uk/) and [JT McCrone](https://jtmccr1.github.io/).
It uses lineages from [Rambaut et al.](https://www.biorxiv.org/content/10.1101/2020.04.17.046086v1).


### Citing ``pangolin``

There is a publication in prep for ``pangolin``, but in the meantime please to link to this github [github.com/hCoV-2019/pangolin](github.com/hCoV-2019/pangolin) if you have used ``pangolin`` in your research. 

### References

The following external software is run as part of pangolin:

[iqtree](http://www.iqtree.org/#download)

L.-T. Nguyen, H.A. Schmidt, A. von Haeseler, B.Q. Minh (2015) IQ-TREE: A fast and effective stochastic algorithm for estimating maximum likelihood phylogenies.. Mol. Biol. Evol., 32:268-274. https://doi.org/10.1093/molbev/msu300

D.T. Hoang, O. Chernomor, A. von Haeseler, B.Q. Minh, L.S. Vinh (2018) UFBoot2: Improving the ultrafast bootstrap approximation. Mol. Biol. Evol., 35:518–522. https://doi.org/10.1093/molbev/msx281

Stéphane Guindon, Jean-François Dufayard, Vincent Lefort, Maria Anisimova, Wim Hordijk, Olivier Gascuel, New Algorithms and Methods to Estimate Maximum-Likelihood Phylogenies: Assessing the Performance of PhyML 3.0, Systematic Biology, Volume 59, Issue 3, May 2010, Pages 307–321, https://doi.org/10.1093/sysbio/syq010

[mafft](https://mafft.cbrc.jp/alignment/software/)

Katoh, Standley 2013 (Molecular Biology and Evolution 30:772-780)
MAFFT multiple sequence alignment software version 7: improvements in performance and usability.
(outlines version 7)

[snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.
