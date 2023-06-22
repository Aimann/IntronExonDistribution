## Path to your gtf file
gtf_file = "gencode.v41.basic.annotation.gtf"

## Path to your bam files
path_to_bams = "bams/{condition}.bam"

## Transcript types; can do multiple types by separating with a space
genetype = "protein_coding lncRNA"

## Get all the samples
samples = sorted(glob_wildcards(path_to_bams).condition)

rule all:
    input:
        "intron_exon_rpkm.txt",
        "intron_exon_rpkm_summary.txt",
	    "intron_exon_rpkm_reformat.txt"
## Get the exon and intron coordinates and generates a gtf file
rule get_exon_intron_coords:
    input:
        gtf=gtf_file
    output:
        "introns_exons.gtf",
	    "gene_types.txt"
    shell:
        "python scripts/intronExon.py --gtf {input.gtf} --genetype {genetype}"
## Count the reads in the introns and the exons
rule count_overlapping_exons:
    input:
        bams=expand(path_to_bams,condition  = samples),
        exons="introns_exons.gtf"
    output:
        "intron_exon_counts.txt"
    shell:
        "featureCounts -s 2 -O --largestOverlap -t exon --fracOverlap 0.2 -a {input.exons} -o {output} {input.bams}"
## Convert the counts to RPKM values
rule counts_to_rpkm:
    input:
        count_file="intron_exon_counts.txt",
	gene_type_file="gene_types.txt"
    output:
        "intron_exon_rpkm.txt",
        "intron_exon_rpkm_summary.txt",
	    "intron_exon_counts_reformat.txt",
	    "intron_exon_rpkm_reformat.txt"
    shell:
        "python scripts/calculateRPKM.py --counts {input.count_file} --genetypes {input.gene_type_file}"
