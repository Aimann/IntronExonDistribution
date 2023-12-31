# IntronExonDistribution


## About

Small set of python scripts to map reads from mRNA sequencing data to the introns and exons of the longest protein coding gene transcripts in a given GTF file

## Installation

Dependencies can be installed using conda:

```bash
conda env create -f intronexon.yaml
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
```bash
## Path to your gtf file
gtf_file = "/path/to/gtf/gencode.v41.basic.annotation.gtf"
```
2. Update the Snakefile with the path to your bam files. Can be used to grab all bam files in a directory.
```bash
## Path to your bam files
path_to_bams = "/path/to/bams/{condition}.bam"
```
3. Navigate to the directory with your Snakefile and run it with:
```bash
snakemake --cores 1
```

### Output files
* `intron_exon_counts.txt` contains the intronic and exonic read counts for each gene in the provided GTF file
* `intron_exon_counts_reformat.txt` contains the intronic and exonic read counts for each gene in the provided GTF file formatted in a gene x exon...intron table
* `intron_exon_rpkm.txt` contains the intronic and exonic RPKM values for each gene in the provided GTF file
* `intron_exon_rpkm_reformat.txt` contains the intronic and exonic RPKM values for each gene in the provided GTF file formatted in a gene x exon...intron table. Also contains the caluclated ratio of introns/exons for each gene. For genes with no introns it labels these as "no intron"
* `intron_exon_rpkm_summary.txt` contains a summary of the intronic and exonic mapping reads for each feature type from each input bam file.
* `introns_exons.gtf` contains the gene coordinates for the introns and exons for each of the longest transcripts in the provided GTF file.
* `gene_types.txt` contains the gene types corresponding to each gene from the input gtf


