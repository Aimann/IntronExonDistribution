#!/usr/bin/env python3
#Aidan Manning 6.14.23
#Script to determine sites of tRNA modficiation from sequencing data and generates some plots
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import csv

class exonIntron(object):
    
    def extract_coordinates(self, gtf_file):
        gene_transcripts = {}
        with open(gtf_file, 'r') as file:
            for line in file:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if parts[2] == 'exon':
                    attributes = dict(item.strip().split(' ') for item in parts[8].strip().split(';') if item.strip())
                    if '_PAR_' in attributes['gene_id'].strip('"'):
                        pass
                    elif '_PAR_' not in attributes['transcript_id'].strip('"') and attributes['gene_type'].strip('"') == args.genetype:
                        gene_id = attributes['gene_name'].strip('"')
                        transcript_id = attributes['transcript_name'].strip('"')
                        chrom = parts[0]
                        strand = parts[6]
                        exon_start = int(parts[3])
                        exon_end = int(parts[4])
                        if gene_id not in gene_transcripts:
                            gene_transcripts[gene_id] = {}
                        if transcript_id not in gene_transcripts[gene_id]:
                            gene_transcripts[gene_id][transcript_id] = []
                        gene_transcripts[gene_id][transcript_id].append((chrom, exon_start, exon_end, strand))
                    else:
                        pass

        return gene_transcripts

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)

    ap.add_argument('--gtf', required=True, help='gtf file containing gene features')
    ap.add_argument('--genetype', required=True, help='type of genes to extract from gtf file; default: protein_coding', default='protein_coding')
    args = ap.parse_args()
    
    ## Generates exonIntron object
    exonintron = exonIntron()

    startTime = datetime.now()
    ## Reads in gff file and generates a database
    print('Reading in gtf file...')

    genes = exonintron.extract_coordinates(args.gtf)
    genedata = []

    for gene_id, transcripts in genes.items():
            
        longest_transcript = max(transcripts, key=lambda x: sum(y[2] - y[1] for y in transcripts[x]))
        exons = transcripts[longest_transcript]
        if exons[0][3] == '+':
            introns = [(exons[i][0] ,exons[i][2] + 1, exons[i + 1][1] - 1, exons[i][3]) for i in range(len(exons) - 1)]
        elif exons[0][3] == '-':
            introns = [(exons[i][0] ,max(exons[i + 1][1] - 1, 0), max(exons[i][2] + 1, 0), exons[i][3]) for i in range(len(exons) - 1)]

        for exon_num, (chrom, exon_start, exon_end, strand) in enumerate(exons, start=1):
            infocol = ''.join([f'gene_id "{gene_id}_exons"; transcript_id "{longest_transcript}_exons"; exon_number "{exon_num}";'])
            genedata.append([chrom, '.', 'exon', exon_start, exon_end, '.', strand, '.', infocol])

        for intron_num, (chrom, intron_start, intron_end, strand) in enumerate(introns, start=1):
            infocolintron = ''.join([f'gene_id "{gene_id}_introns"; transcript_id "{longest_transcript}_introns"; intron_number "{intron_num}";'])
            genedata.append([chrom, '.', 'exon', intron_start, intron_end, '.', strand, '.', infocolintron])
            # introndata.append([chrom, '.', 'intron', intron_start, intron_end, '.', strand, '.', infocolintron])

    genedf = pd.DataFrame(genedata, columns=['chrom', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    genedf.to_csv('introns_exons.gtf', index=False, sep='\t', header=False,quoting=csv.QUOTE_NONE, quotechar="",  escapechar="\\")

    print('\nTime elasped: ', datetime.now() - startTime)