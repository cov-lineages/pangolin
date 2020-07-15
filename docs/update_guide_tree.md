# Updating the guide tree and alignment with new data

### Requirements

1. A csv with taxon names and lineage

2. A fasta file (doesn't need to be an alignment but can be) with the corresponding sequences.

3. A metadata csv from ``grapevine``

### Usage

1. Activate the environment ``conda activate pangolin``
2. Run the following:

```
snakemake 
        --snakefile prepare_package_data.smk    \
        --config                                \
                metadata=your_new_metadata.csv  \
                fasta=your_new_fasta.fasta      \
				global_tree=global_tree.nexus   \
                outdir=where/to/put/your/data   \
                lineages=lineages.csv \
        --cores 2
```
