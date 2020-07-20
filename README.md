# pangolin

**P**hylogenetic **A**ssignment of **N**amed **G**lobal **O**utbreak **LIN**eages

<img src="https://github.com/cov-lineages/pangolin/blob/master/docs/logo.png" width="300">


>pangolin now comes in web-application form thanks to the [Centre for Genomic Pathogen and Surveillance](https://www.pathogensurveillance.net/)! Find it here at https://pangolin.cog-uk.io/.


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
  * [Citing <strong>pangolin</strong>](#citing-pangolin)
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

> Note: we recommend using pangolin in the conda environment specified in the ``environment.yml`` file as per the instructions above. If you can't use conda for some reason, bear in mind the data files are now hosted in a separate repository at [``cov-lineages/lineages``](https://github.com/cov-lineages/lineages.git) and you will need to pip install that alongside the other dependencies for pangolin (details found in [``environment.yml``](https://github.com/cov-lineages/pangolin/blob/master/environment.yml)).

### Check the install worked

Type (in the pangolin environment):

```
pangolin -v
pangolin -pv
pangolin -lv
```
and you should see the versions of pangolin , pangoLEARN and lineages data release printed respectively.

### Updating pangolin

> Note: Even if you have previously installed <strong>pangolin</strong>, as it is being worked on intensively, we recommend you check for updates before running.

To update:

1. ``conda activate pangolin``
2. ``git pull`` \
pulls the latest changes from github
3. ``python setup.py install`` \
re-installs pangolin.
4. ``conda env update -f environment.yml`` \
updates the conda environment (you're unlikely to need to do this, but just in case!)
5. ``pip install git+https://github.com/cov-lineages/pangoLEARN.git --upgrade`` \
updates if there is a new data release
6. ``pip install git+https://github.com/cov-lineages/lineages.git --upgrade`` \
updates if there is a new data release, this is the legacy data repo and is unlikely to have tagged releases in the future

### Usage

1. Activate the environment ``conda activate pangolin``
2. Run ``pangolin <query>``, where ``<query>`` is the name of your input file.

```
pangolin: Phylogenetic Assignment of Named Global Outbreak LINeages

positional arguments:
  query                 Query fasta file of sequences to analyse.

optional arguments:
  -h, --help            show this help message and exit
  -o OUTDIR,
  --outdir OUTDIR
                        Output directory. Default: current working directory
  --outfile OUTFILE     Optional output file name. Default: lineage_report.csv
  -d DATA, --data DATA  Data directory minimally containing a fasta alignment
                        and guide tree
  -n, --dry-run         Go through the motions but don't actually run
  --tempdir TEMPDIR     Specify where you want the temp stuff to go. Default:
                        $TMPDIR
  --no-temp             Output all intermediate files, for dev purposes.
  --max-ambig MAXAMBIG  Maximum proportion of Ns allowed for pangolin to
                        attempt assignment. Default: 0.5
  --min-length MINLEN   Minimum query length allowed for pangolin to attempt
                        assignment. Default: 10000
  --panGUIlin           Run web-app version of pangolin
  --assign-using-tree   LEGACY: Use original phylogenetic assignment methods
                        with guide tree. Note, will be significantly slower
                        than pangoLEARN
  --write-tree          Output a phylogeny for each query sequence placed in
                        the guide tree. Only works in combination with legacy
                        `--assign-using-tree`
  -t THREADS,
  --threads THREADS     Number of threads
  -p, --include-putative
                        Include the bleeding edge lineage definitions in
                        assignment
  --verbose             Print lots of stuff to screen
  -v, --version         show program's version number and exit
  -lv, --lineages-version
                        show lineages's version number and exit
  -pv, --pangoLEARN-version
                        show pangoLEARN's version number and exit
  ```

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

Legacy phylogenetics output example:

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

Given cov-lineages is relatively slow evolving for an RNA virus and there is still not a huge amount of diversity, missing or ambiguous data at key residues may lead to incorrect placement within the guide tree. We have a filter in place that by default with not call a lineage for any sequence with >50% N-content, but this can be made more conservative with the command line option `--max-ambig`.

### Source data

<strong>pangolin</strong> runs a multinomial logistic regression model trained against lineage assignments based on GISAID data.

Legacy <strong>pangolin</strong> runs using a guide tree and alignment hosted at [``cov-lineages/lineages``](https://github.com/cov-lineages/lineages.git). Some of this data is sourced from GISAID, but anonymised and encrypted to fit with guidelines. Appropriate permissions have been given and acknowledgements for the teams that have worked to provide the original SARS-CoV-2 genome sequences to GISAID are also hosted in [``cov-lineages.org``](https://cov-lineages.org/gisaid_acknowledgements.html).


### Authors

Pangolin was created by [Áine O'Toole](https://aineotoole.co.uk/), [JT McCrone](https://jtmccr1.github.io/) and [Emily Scher](http://www.scherscherscher.com/).
It uses lineages from [Rambaut et al. 2020](https://www.nature.com/articles/s41564-020-0770-5).


### Citing <strong>pangolin</strong>

There is a publication in prep for <strong>pangolin</strong>, but in the meantime please to link to this github [github.com/cov-lineages/pangolin](github.com/cov-lineages/pangolin) if you have used <strong>pangolin</strong> in your research. 

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

[minimap2](https://github.com/lh3/minimap2)

Heng Li, Minimap2: pairwise alignment for nucleotide sequences, Bioinformatics, Volume 34, Issue 18, 15 September 2018, Pages 3094–3100, https://doi.org/10.1093/bioinformatics/bty191

[snakemake](https://snakemake.readthedocs.io/en/stable/index.html)

Köster, Johannes and Rahmann, Sven. “Snakemake - A scalable bioinformatics workflow engine”. Bioinformatics 2012.
