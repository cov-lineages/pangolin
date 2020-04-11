# Updating the guide tree and alignment with new data

### Requirements

1. A csv with taxon names, lineage and whether or not it's a representative sequence.

Example:

```
name,lineage,bootstrap,amiguity,representative
Hangzhou/ZJU-08/2020|EPI_ISL_416473||China|Hangzhou||2020-01-26,A,99/67,0,1
Wuhan/HBCDC-HB-02/2020|EPI_ISL_412978||China|Hubei|Wuhan|2020-01-17,A,96,0,1
Wuhan/HBCDC-HB-04/2020|EPI_ISL_412980||China|Hubei|Wuhan|2020-01-18,A,None,0,1
```
2. A fasta file (doesn't need to be an alignment but can be) with the corresponding sequences.

### Usage

1. Activate the environment ``conda activate pangolin-env``
2. Run the following:

```
snakemake 
        --snakefile prepare_package_data.smk    \
        --config                                \
                metadata=your_new_metadata.csv  \
                fasta=your_new_fasta.fasta      \
                outdir=where/to/put/your/data   \
        --cores 2
```
