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

Grab the most recent gencode gtf file

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_41/gencode.v41.basic.annotation.gtf.gz
```

unzip this file
```bash
gunzip gencode.v41.basic.annotation.gtf.gz
```

## Usage

This script extracts the longest isoform of each gene in a given gtf file tagged as 'protein_coding' (default) and splits it into its respective introns and exons. It then determines the number of reads mapping to each of these transcript features only being assigned to the feature it overlaps the most. 

### Activating the environment

```bash
conda activate intronexon
```

### Getting the exonic and intronic mapping reads from a set of RNA-seq bam files

1. Update the Snakefile with the full path to your Gencode GTF file.
2. Update the Snakefile with the path to your bam files. Can be used to grab all bam files in a directory.
3. Navigate to the directory with your Snakefile and run it with:

```bash
snakemake
```
### Output files
* `intron_exon_counts.tx` contains the intronic and exonic read counts for each gene in the provided GTF file
* `intron_exon_rpkm.txt` contains the intronic and exonic RPKM values for each gene in the provided GTF file
* `intron_exon_rpkm_summary.txt` contains a summary of the intronic and exonic mapping reads from each input bam file.
* `introns_exons.gtf` contains the intronic and exonic RPKM values for each gene in the provided GTF file