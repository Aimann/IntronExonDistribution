# IntronExonDistribution


## About

Small set of python scripts to map reads from mRNA sequencing data to the introns and exons of the longest protein coding gene transcripts in a given GTF file

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f intronexon.yml
```

## Input files

### Gencode GTF File

Can just grab the most recent gencode gtf file using the following

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gtf.gz
```

and then unzip this file
```bash
gunzip gencode.v41.basic.annotation.gtf.gz
```

## Usage

### Activating the environment

```bash
conda activate intronexon
```

### Getting the exonic and intronic mapping reads from a set of RNA-seq bam files

This script extracts the longest isoform of each gene in a given gtf file tagged as 'protein_coding' (default) and splits it into its respective introns and exons. It then determines the number of reads mapping to each of these transcript features only being assigned to the feature it overlaps the most. 

```bash
snakemake
```

* `--cov` is the path to the tRAX coverage file you want to analyse modifications from.
* `--alignments` is the path to the stockholm file containing the mature tRNA alignments from tRAX makedb.py; dbname-trnaalign.stk.
* `--o`  is the prefix for the output files.
* (optional flags)
* `--minreads` is the minumum number of read coverage of a base required to go into analysis; default=20.
* `--org` organism of interest (euk, bact, arch, mito); default='euk'.
* `--plot` whether or not to generate a bar plot showing the Sprinzl positions containing instances of isodecoder variation.