#!/usr/bin/env python3
#Aidan Manning 6.14.23
#Script to calculate the RPKM values of the intron and exon counts
import argparse
import sys
import pandas as pd
import os
from datetime import datetime

class exonIntronRPKM(object):
    
    def calculateRPKM(self, counts_file):
        counts = pd.read_csv(counts_file, sep='\t', header=0, index_col=0, skiprows=[0])
        counts = counts.fillna(0)
        to_drop = ['Chr', 'Start', 'End', 'Length', 'Strand']
        for col in list(counts.columns)[5:]:
            counts[col + '|RPKM'] = round(counts[col] / ((counts['Length'] / 1000) * (counts[col].sum() / 1000000)), 2).fillna(0)
            to_drop.append(col)
        rpkmcounts = counts.drop(to_drop, axis=1)

        outdata = []
        for col in rpkmcounts.columns:
            introncounts = rpkmcounts.filter(like='intron', axis=0)
            exoncounts = rpkmcounts.filter(like='exon', axis=0)
            totalcounts = introncounts[col].sum() + exoncounts[col].sum()
            outdata.append([col, introncounts[col].sum(), exoncounts[col].sum(), round((introncounts[col].sum() / totalcounts)*100), round((exoncounts[col].sum() / totalcounts)*100)])
            # print(col + '\t' + str(introncounts[col].sum()) + '\t' + str(exoncounts[col].sum()) + '\t' + str(round((introncounts[col].sum() / totalcounts)*100)) + '\t' + str(round((exoncounts[col].sum() / totalcounts)*100)) + '\n')
        outdf = pd.DataFrame(outdata, columns=['Sample', 'Total Introns RPKM', 'Total Exons RPKM', '% Introns', '% Exons'])

        return rpkmcounts, outdf

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)

    ap.add_argument('--counts', required=True, help='counts of introns and exons')
    args = ap.parse_args()
    
    ## Generates exonIntron object
    rpkm = exonIntronRPKM()

    startTime = datetime.now()
    ## Reads in gff file and generates a database
    print('Reading calculating RPKM values...')

    rpkm_values, out_data = rpkm.calculateRPKM(args.counts)

    rpkm_values.to_csv('intron_exon_rpkm.txt', sep='\t', index=True)
    out_data.to_csv('intron_exon_rpkm_summary.txt', sep='\t', index=False)

    print('\nTime elasped: ', datetime.now() - startTime)