#!/usr/bin/env python3
import pandas as pd
import argparse
import sys
from pybedtools import BedTool
import gzip
import numpy as np


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--pos_selected',
                        help='A file with columns chrom, start, end, category. The file should be tab delimited. '
                             'Coordinates must be in hg19')
    parser.add_argument('--neg_selected',
                        help='A file with columns chrom, start, end, category. The file should be tab delimited. '
                             'Coordinates must be in hg19')
    parser.add_argument('--deserts',
                        help='A file with columns chrom, start, end, category. The file should be tab delimited. '
                             'Coordinates must be in hg19')
    parser.add_argument('--baseline', help='LDSC baseline annotation files')
    parser.add_argument('--annot-file', type=str, help='the name of the annot file to output.')
    args = parser.parse_args()
    # load annotations
    pos_selected = pd.read_csv(args.pos_selected, sep='\t', names=['chrom', 'start', 'end', 'category'])
    pos_selected.start = np.where(pos_selected.start.values < 0, 0, pos_selected.start.values)
    pos_selected.sort_values(['chrom', 'start'], inplace=True)

    neg_selected = pd.read_csv(args.neg_selected, sep='\t', names=['chrom', 'start', 'end', 'category'])
    neg_selected.start = np.where(neg_selected.start.values < 0, 0, neg_selected.start.values)
    neg_selected.sort_values(['chrom', 'start'], inplace=True)

    deserts = pd.read_csv(args.deserts, sep='\t', names=['chrom', 'start', 'end', 'category'])
    deserts.start = np.where(deserts.start.values < 0, 0, deserts.start.values)
    deserts.sort_values(['chrom', 'start'], inplace=True)
    # load baseline model
    baseline = pd.read_csv(args.baseline, sep='\t', header=0)
    categories = baseline.columns[5:]
    for cat in categories:
        baseline.drop(cat, axis=1, inplace=True)
    baseline.BP = baseline.BP.astype(int)
    # get bed file
    baselinebed = BedTool([[chrom, bp - 1, bp] for chrom, bp in baseline.loc[:, ['CHR', "BP"]].values])
    # get current annotation
    category_pos_selected = pos_selected.category.unique()[0]
    category_neg_selected = neg_selected.category.unique()[0]
    category_deserts = deserts.category.unique()[0]

    pos_selected.drop('category', axis=1, inplace=True)
    neg_selected.drop('category', axis=1, inplace=True)
    deserts.drop('category', axis=1, inplace=True)
    # calculate overlap between annotation and baseline model
    pos_selected_bed = BedTool.from_dataframe(pos_selected)
    neg_selected_bed = BedTool.from_dataframe(neg_selected)
    deserts_bed = BedTool.from_dataframe(deserts)

    overlap_pos_selected = baselinebed.intersect(pos_selected_bed, c=True).to_dataframe()
    overlap_neg_selected = baselinebed.intersect(neg_selected_bed, c=True).to_dataframe()
    overlap_deserts = baselinebed.intersect(deserts_bed, c=True).to_dataframe()

    # modify baseline annotations based on overlap
    baseline[category_pos_selected] = overlap_pos_selected.name.values
    baseline[category_neg_selected] = overlap_neg_selected.name.values
    baseline[category_deserts] = overlap_deserts.name.values

    # categories = baseline.columns[5:-1]
    # for cat in categories:
    #     #vals = baseline.loc[:, cat].values
    #     #vals = np.where(overlap.name.values == 1, vals, 0)
    #     #baseline[f'{cat}_{category}'] = vals
    #     baseline.drop(cat, axis=1, inplace=True)
    baseline.BP = baseline.BP.astype(int)
    # save
    if args.annot_file.endswith('.gz'):
        with gzip.open(args.annot_file, 'wb') as f:
            baseline.to_csv(f, sep="\t", index=False)
    else:
        baseline.to_csv(args.annot_file, sep="\t", index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
