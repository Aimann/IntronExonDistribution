#!/usr/bin/env python3
#Aidan Manning 6.14.23
#Script to calculate the RPKM values of the intron and exon counts
import argparse
import sys
import pandas as pd
import os
from datetime import datetime
import warnings
warnings.filterwarnings("ignore")

class exonIntronRPKM(object):

    def __init__(self, genetypes):
        
        self.genetypeD = dict(zip(genetypes['index'], genetypes['gene_type']))
        self.genetypes = genetypes['gene_type'].drop_duplicates().tolist()

        return

    def reformatCounts(self, cts):

        raw_counts = pd.read_csv(cts, sep='\t', header=0, index_col=0, skiprows=[0])
        raw_counts = raw_counts.drop(['Chr', 'Start', 'End', 'Length', 'Strand'], axis=1)
        raw_counts_exons = raw_counts.filter(like='exon', axis=0).reset_index()
        raw_counts_introns = raw_counts.filter(like='intron', axis=0).reset_index()
        raw_counts_exons['Geneid'] = raw_counts_exons['Geneid'].str.replace('_exons', '')
        raw_counts_exons = raw_counts_exons.set_index('Geneid').rename(columns={samplename: samplename + '_exons' for samplename in raw_counts_introns.columns if samplename != 'Geneid'})
        raw_counts_introns['Geneid'] = raw_counts_introns['Geneid'].str.replace('_introns', '')
        raw_counts_introns = raw_counts_introns.set_index('Geneid').rename(columns={samplename: samplename + '_introns' for samplename in raw_counts_introns.columns if samplename != 'Geneid'})
        reformatted_raw_counts = pd.concat([raw_counts_exons, raw_counts_introns], axis=1, sort=True).fillna('no introns')
        colorder = reformatted_raw_counts.columns.tolist()
        reformatted_raw_counts['gene_type'] = reformatted_raw_counts.index.map(self.genetypeD)
        reformatted_raw_counts = reformatted_raw_counts[['gene_type'] + colorder]

        return reformatted_raw_counts

    def calculateRPKM(self, counts_file):
        counts = pd.read_csv(counts_file, sep='\t', header=0, index_col=0, skiprows=[0])
        counts = counts.fillna(0)
        to_drop = ['Chr', 'Start', 'End', 'Length', 'Strand']
        for col in list(counts.columns)[5:]:
            counts[col + '|RPKM'] = round(counts[col] / ((counts['Length'] / 1000) * (counts[col].sum() / 1000000)), 2).fillna(0)
            to_drop.append(col)
        rpkmcounts = counts.drop(to_drop, axis=1)
        colorder = rpkmcounts.columns.tolist()
        rpkmcounts['ID'] = rpkmcounts.index.str.replace('_exons', '').str.replace('_introns', '')
        rpkmcounts['gene_type'] = rpkmcounts['ID'].map(self.genetypeD)
        rpkmcounts = rpkmcounts.drop(['ID'], axis=1)

        outdata = []
        for col in rpkmcounts.columns[:-1]:
            for genetype in self.genetypes:
                rpkmcounts_genetype = rpkmcounts[rpkmcounts['gene_type'] == genetype]
                introncounts = rpkmcounts_genetype.filter(like='intron', axis=0)
                exoncounts = rpkmcounts_genetype.filter(like='exon', axis=0)
                totalcounts = introncounts[col].sum() + exoncounts[col].sum()
                outdata.append([col, introncounts[col].sum(), exoncounts[col].sum(), round((introncounts[col].sum() / totalcounts)*100, 2), round((exoncounts[col].sum() / totalcounts)*100, 2), genetype])
            # print(col + '\t' + str(introncounts[col].sum()) + '\t' + str(exoncounts[col].sum()) + '\t' + str(round((introncounts[col].sum() / totalcounts)*100)) + '\t' + str(round((exoncounts[col].sum() / totalcounts)*100)) + '\n')
        outdf = pd.DataFrame(outdata, columns=['Sample', 'Total Introns RPKM', 'Total Exons RPKM', '% Introns', '% Exons', 'Gene Class'])

        rpkmcounts = rpkmcounts[['gene_type'] + colorder]

        return rpkmcounts, outdf

    def reformatRPKM(self, rpkm):

        rpkm_nogene = rpkm[rpkm.columns[1:]]

        samplelist = rpkm_nogene.columns.tolist()
        rpkm_exons = rpkm_nogene.filter(like='exon', axis=0).reset_index()
        rpkm_introns = rpkm_nogene.filter(like='intron', axis=0).reset_index()
        rpkm_exons['Geneid'] = rpkm_exons['Geneid'].str.replace('_exons', '')
        rpkm_exons = rpkm_exons.set_index('Geneid').rename(columns={samplename: samplename + '_exons' for samplename in rpkm_introns.columns if samplename != 'Geneid'})
        rpkm_introns['Geneid'] = rpkm_introns['Geneid'].str.replace('_introns', '')
        rpkm_introns = rpkm_introns.set_index('Geneid').rename(columns={samplename: samplename + '_introns' for samplename in rpkm_introns.columns if samplename != 'Geneid'})
        reformatted_rpkm = pd.concat([rpkm_exons, rpkm_introns], axis=1, sort=True).fillna('no introns')
        colorder = reformatted_rpkm.columns.tolist()
        reformatted_rpkm['gene_type'] = reformatted_rpkm.index.map(self.genetypeD)

        reformatted_rpkm = reformatted_rpkm[['gene_type'] + colorder]
        no_introns = reformatted_rpkm[reformatted_rpkm[colorder[-1]] == 'no introns']
        reformatted_rpkm_introncontain = reformatted_rpkm[reformatted_rpkm[colorder[-1]] != 'no introns']

        for sample in samplelist:
            reformatted_rpkm_introncontain[sample + '|Intron_Exon_Ratio'] = round(reformatted_rpkm_introncontain[sample + '_introns'].astype(float) / reformatted_rpkm_introncontain[sample + '_exons'].astype(float), 2).fillna(0)
            colorder.append(sample + '|Intron_Exon_Ratio')
        rpkm_out = pd.concat([reformatted_rpkm_introncontain, no_introns], axis=0, sort=True).fillna('no introns')
        rpkm_out = rpkm_out[['gene_type'] + colorder]

        return rpkm_out

if __name__ == '__main__':

    ap = argparse.ArgumentParser(prog=os.path.basename(sys.argv[0]), usage=__doc__)

    ap.add_argument('--counts', required=True, help='counts of introns and exons')
    ap.add_argument('--genetypes', required=True, help='genetypes file')
    args = ap.parse_args()
    
    ## Generates exonIntron object
    rpkm = exonIntronRPKM(pd.read_csv(args.genetypes, sep='\t', header=0, index_col=0).reset_index())

    startTime = datetime.now()
    ## Reads in gff file and generates a database
    print('Calculating RPKM values...')

    reformatted_raw_counts = rpkm.reformatCounts(args.counts)
    reformatted_raw_counts.to_csv('intron_exon_counts_reformat.txt', sep='\t', index=True)

    rpkm_values, out_data = rpkm.calculateRPKM(args.counts)

    reformatted_rpkm = rpkm.reformatRPKM(rpkm_values)

    reformatted_rpkm.to_csv('intron_exon_rpkm_reformat.txt', sep='\t', index=True)

    rpkm_values.to_csv('intron_exon_rpkm.txt', sep='\t', index=True)
    out_data.to_csv('intron_exon_rpkm_summary.txt', sep='\t', index=False)

    print('\nTime elasped: ', datetime.now() - startTime)