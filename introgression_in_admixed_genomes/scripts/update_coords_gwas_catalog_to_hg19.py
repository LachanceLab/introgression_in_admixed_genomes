#!/usr/bin/env python3
import pandas as pd
import re
import numpy as np
import sys
import argparse


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--gwas', help='GWAS catalog in .tab format with coordinates referenced against hg38')
    parser.add_argument('--dbsnp_hg19', help='dbSNP database in BED format with coordinates referenced against hg19')
    parser.add_argument('-o', '--output', help='Output file, GWAS catalog with hg19 coordinates')
    args = parser.parse_args()
    gwas = pd.read_csv(args.gwas, sep='\t', header=0)
    columns = ['chrom', 'start', 'end']
    columns.extend(gwas.columns.values[3:].tolist())
    dbsnp = pd.read_csv(args.dbsnp_hg19, sep='\t', names=['chrom', 'start', 'end', 'rsID'], usecols=[0, 1, 2, 3])
    merged = gwas.set_index('SNPS').join(dbsnp.set_index('rsID'))
    merged = merged[~np.isnan(merged.start.values)]
    merged.reset_index(inplace=True)
    merged = merged.loc[:, columns]
    # get autosomes
    merged = merged.loc[[True if re.match('chr[1-9][0-9]?$', chrom) else False for chrom in merged.chrom]]
    # format chromosomes
    merged.chrom = [int(chrom.replace('chr', '')) for chrom in merged.chrom]
    merged.chrom = merged.chrom.astype(int)
    merged.start = merged.start.astype(int)
    merged.end = merged.end.astype(int)
    merged.sort_values(['chrom', 'start', 'end'], inplace=True)
    merged.to_csv(args.output, header=True, sep='\t', index=False)

if __name__ == '__main__':
    main(sys.argv[1:])
