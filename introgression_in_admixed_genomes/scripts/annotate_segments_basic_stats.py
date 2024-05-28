#!/usr/bin/env python3
import sys
import argparse
import pandas as pd
import numpy as np
from pybedtools import BedTool


def calculate_unique_overlap(df):
    """
    Calculate overlap with, e.g., functional elements
    :param df: pd.DataFrame, BedTool intersect dataframe
    :rtype: int, coverage
    """
    if np.all(df.overlap.values == 0):
        return 0
    else:
        segment_length = (df.end - df.start).values[0]
        df1 = BedTool.from_dataframe(df.loc[:, ['chrom', 'start', 'end']].drop_duplicates())
        df2 = BedTool.from_dataframe(df.loc[:, ['chrom_b', 'start_b', 'end_b']].drop_duplicates())
        df = df1.intersect(df2).sort()
        df = df.merge().to_dataframe()
        return sum(df.end - df.start) / segment_length


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--intervals', help='Bed file with introgressed segments to annotated')
    parser.add_argument('-c', '--config', help='Config file for annotations')
    parser.add_argument('--quantitative', action='store_true', help='Of categorical data get only quantitative '
                                                                      'summary stats, e.g., overlap with genes '
                                                                      'but not gene names', default=False)
    parser.add_argument('-o', '--output', help="Output file name")
    args = parser.parse_args()
    # weight mean by overlap
    weighted_mean = lambda b, w: sum(np.array(b) * np.array(w)) / sum(w)
    # weight standard deviation by overlap
    weighted_stdev = lambda b, w: np.sqrt(sum(np.array(w) * (weighted_mean(b, w) - np.array(b)) ** 2) / sum(w))
    # calculate overlap
    median = lambda b, w: np.sort(b)[np.cumsum(np.array(w)[np.argsort(b)])[
        np.cumsum(np.array(w)[np.argsort(b)]) <= np.cumsum(np.array(w)[np.argsort(b)])[-1] // 2].shape[0]]

    df = pd.read_csv(args.intervals, header=0, sep='\t')
    annotations = []
    intervals = BedTool(args.intervals)
    with open(args.config, 'r') as conf:
        # iterate over all annotation files
        for line in conf:
            file, cols, target_col, data_type = line.strip().split('\t')
            columns = df.columns.values.tolist()

            columns.extend(cols.split(','))
            columns.append('overlap')
            # intersect
            anno = intervals.intersect(file, wao=True).to_dataframe(names=columns)

            # numerical annotation --> can calculate mean, median, and standard deviation
            if data_type == 'num':
                anno = anno.groupby(['chrom', 'start', 'end']).agg(target=(target_col, list),
                                                                   weights=('overlap', list))
                avg = []
                stdev = []
                median_vals = []
                distinct = []
                # calculate stats
                for target, weights in zip(anno.target.values, anno.weights.values):
                    if target == [-1]:
                        avg.append(np.nan)
                        stdev.append(np.nan)
                        median_vals.append(np.nan)
                        distinct.append('.')
                    else:
                        target = np.array([t for t in target if t != '.']).astype(float)
                        if target.shape[0] == 0:
                            avg.append(np.nan)
                            stdev.append(np.nan)
                            median_vals.append(np.nan)
                            distinct.append('.')
                            continue
                        avg.append(weighted_mean(target, weights))
                        stdev.append(weighted_stdev(target, weights))
                        median_vals.append(median(target, weights))
                        if len(target) > 250 and args.quantitative:
                            target = np.random.choice(target, 250, replace=False)
                        distinct.append(",".join([str(val) for val in target]))
                anno[f'{target_col}_avg'] = avg
                anno[f'{target_col}_stdev'] = stdev
                anno[f'{target_col}_median'] = median_vals
                anno[f'{target_col}_distinct'] = ['' if isinstance(val, float) and np.isnan(val) else str(val)
                                                  for val in distinct]

                anno.drop(['target', 'weights'], axis=1, inplace=True)

            # categorical annotation such as overlap with genes or regulatory regions --> calculate coverage
            elif data_type == 'str':
                coverage = anno.drop_duplicates(['chrom', 'start', 'end']).apply(lambda row: calculate_unique_overlap(anno[(anno.chrom == row['chrom']) &
                                                                                (anno.start == row['start']) &
                                                                                (anno.end == row['end'])]), axis=1)
                # this removes alternative splice isoforms
                no_overlap = anno[anno.overlap == 0]
                unique_anno = anno[anno['overlap'] != 0]
                # unique_anno.drop_duplicates(['chrom_b', 'start_b', 'end_b'], inplace=True)
                anno = pd.concat([unique_anno, no_overlap]).sort_values(['chrom', 'start', 'end'])
                anno = anno.groupby(['chrom', 'start', 'end']).agg(target=(target_col, list))
                anno[f'{target_col}_coverage'] = coverage.values
                if args.quantitative:
                    anno.drop(columns='target', inplace=True)
                else:
                    anno.rename({'target': target_col}, axis=1, inplace=True)
            annotations.append(anno)
    conf.close()
    df = df.set_index(['chrom', 'start', 'end']).join(pd.concat(annotations, axis=1))
    df.reset_index(drop=False, inplace=True)
    df.start = df.start.astype(int)
    df.end = df.end.astype(int)
    df.to_csv(args.output, sep='\t', header=True, index=False, na_rep='.')


if __name__ == '__main__':
    main(sys.argv[1:])
