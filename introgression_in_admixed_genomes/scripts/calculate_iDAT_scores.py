#!/usr/bin/env python3
import argparse
import sys
import pandas as pd
import numpy as np
import multiprocessing as mp
from itertools import islice, takewhile, repeat
import re

"""
This script implements the iDAT score from Iman Hamid, Katharine L Korunes, Sandra Beleza, 
Amy Goldberg (2021) Rapid adaptation to malaria facilitated by admixture in the human population of Cabo Verde
 eLife 10:e63177. https://doi.org/10.7554/eLife.63177
"""


def dat_score(df, ref_point, ref_ancestry, dist_step_size=10000, min_cov=0.1):
    """
    Calculate DAT score for given ancestry at specified point
    :param df: pd.DataFrame, ancestry tracks
    :param ref_point: int, genomic position
    :param ref_ancestry: float, current reference ancestry
    :param dist_step_size: int, increment with which to move away from reference point
    :param min_cov: float, minimum DAT score considered for calculating AUC
    :return: np.array-like, np.array-like, distances and corresponding DAT scores
    """
    # get haplotypes overlapping reference point with given ancestry
    haplotypes = df[(df.start <= ref_point) & (df.end > ref_point) & (df.la == ref_ancestry)]
    #  count haplotypes
    n_haplotypes = haplotypes.shape[0]
    if n_haplotypes == 0:
        return np.array([0]), np.array([1e-20])
    dat_scores = [1]
    # increment with which we'll increase the distance
    distances = [0, dist_step_size]
    # calculate overlap to right and left
    overlap = haplotypes[(haplotypes.start <= ref_point - distances[-1]) |
                         (haplotypes.end > ref_point + distances[-1])].shape[0]
    # DAT score
    dat_scores.append((overlap / n_haplotypes) ** 2)
    while dat_scores[-1] >= min_cov:
        # reach start and end of chromosome
        if ref_point - (distances[-1] + dist_step_size) < df.start.min() and\
           ref_point + distances[-1] + dist_step_size > df.end.max():
            break
        # increase distance
        distances.append(distances[-1] + dist_step_size)
        # calculate overlap to right and left

        overlap = haplotypes[(haplotypes.start <= ref_point - distances[-1]) |
                             (haplotypes.end > ref_point + distances[-1])].shape[0]
        # get DAT score
        dat_scores.append((overlap / n_haplotypes) ** 2)
    return np.array(distances) / 1e6, np.array(dat_scores)


def calculate_idat_score(args):
    """
    Calculate iDAT score
    :param args: tuple, pd.DataFrame with ancestry tracks, list-like with reference points,
                        int specifying the step size used to move away from reference point,
                        and float specifying minimal DAT score
    :return: array-like, array-like, reference points and corresponding iDAT scores
    """
    df, points, dist_step_size, min_cov = args
    idat_scores = []
    for point in points:
        # calculate DAT
        distances_afr, ratios_afr = dat_score(df, point, ref_ancestry='AFR', dist_step_size=dist_step_size,
                                              min_cov=min_cov)
        distances_eur, ratios_eur = dat_score(df, point, ref_ancestry='EUR', dist_step_size=dist_step_size,
                                              min_cov=min_cov)
        distances_eas, ratios_eas = dat_score(df, point, ref_ancestry='EAS', dist_step_size=dist_step_size,
                                              min_cov=min_cov)

        # get AUC and add pseudo count
        dat_afr = np.trapz(ratios_afr, distances_afr) + 1e-10
        dat_eur = np.trapz(ratios_eur, distances_eur) + 1e-10
        dat_eas = np.trapz(ratios_eas, distances_eas) + 1e-10
        # calculate iDAT relative to the other two ancestries
        idat_afr = np.log(dat_afr / (dat_eur + dat_eas))
        idat_eur = np.log(dat_eur / (dat_afr + dat_eas))
        idat_eas = np.log(dat_eas / (dat_afr + dat_eur))
        idat_scores.append((idat_afr, idat_eur, idat_eas))
    return points, idat_scores


def get_idat_scores(args):
    local_ancestry = pd.concat([pd.read_csv(flare, sep='\t', header=0) for flare in args.local_ancestry])

    chrom_sizes = pd.read_csv(args.genomefile, sep='\t', names=['chrom', 'size'])
    c_chrom_size = chrom_sizes.loc[chrom_sizes.chrom == f'chr{args.chrom}', 'size'].values[0]
    # get points at which to calculate DAT
    n_sites = np.ceil(args.number_total_sites * c_chrom_size / chrom_sizes['size'].sum()).astype(int)
    points = iter(np.random.randint(0, c_chrom_size, n_sites))
    chunks = takewhile(bool, (list(islice(points, 16)) for _ in repeat(None)))

    ready_to_map = [(local_ancestry, chunk, args.dist_step_size, args.min_cov_dat) for chunk in chunks]
    pool = mp.Pool(processes=args.threads)
    results = pool.map(calculate_idat_score, ready_to_map)
    pool.close()
    pool.join()
    points = np.concatenate([res[0] for res in results])
    idat_scores = np.vstack([res[1] for res in results]).T
    idat_scores_AFR = idat_scores[0][np.argsort(points)]
    idat_scores_EUR = idat_scores[1][np.argsort(points)]
    idat_scores_EAS = idat_scores[2][np.argsort(points)]
    points = np.sort(points)
    with open(args.output, 'w') as out:
        for pos, idat_afr, idat_eur, idat_eas in zip(points, idat_scores_AFR, idat_scores_EUR,idat_scores_EAS):
            out.write(f'chr{args.chrom}\t{pos}\t{idat_afr}\t{idat_eur}\t{idat_eas}\n')
    out.close()


def get_standardized_idat_scores_windows(args):
    df = pd.concat([pd.read_csv(idat, sep='\t', names=['chrom', 'pos', 'idat_AFR', 'idat_EUR', 'idat_EAS'])
                    for idat in args.idat_scores])
    df.sort_values(['chrom', 'pos'], inplace=True)
    chrom_sizes = pd.read_csv(args.genomefile, sep='\t', names=['chrom', 'size'])
    # standardize to zero mean and standard deviation of 1
    df['standardized_idat_AFR'] = (df.idat_AFR.values - df.idat_AFR.mean()) / df.idat_AFR.std()
    df['standardized_idat_EUR'] = (df.idat_EUR.values - df.idat_EUR.mean()) / df.idat_EUR.std()
    df['standardized_idat_EAS'] = (df.idat_EAS.values - df.idat_EAS.mean()) / df.idat_EAS.std()
    outfile = open(args.output, 'w')
    # calculate standardized iDAT for windows and write to file
    for chrom in df.chrom.unique():
        pos = 0
        c_chrom_size = chrom_sizes.loc[chrom_sizes.chrom == chrom, 'size'].values[0]
        points = df.loc[df.chrom == chrom, 'pos'].values
        standardized_idat_AFR = df.loc[df.chrom == chrom, 'standardized_idat_AFR'].values
        standardized_idat_EUR = df.loc[df.chrom == chrom, 'standardized_idat_EUR'].values
        standardized_idat_EAS = df.loc[df.chrom == chrom, 'standardized_idat_EAS'].values
        while pos < c_chrom_size:
            avg_idat_afr = np.mean(standardized_idat_AFR[(points >= pos) & (points < args.windowsize + pos)])
            avg_idat_eur = np.mean(standardized_idat_EUR[(points >= pos) & (points < args.windowsize + pos)])
            avg_idat_eas = np.mean(standardized_idat_EAS[(points >= pos) & (points < args.windowsize + pos)])
            outfile.write(f'{chrom}\t{pos}\t{min([c_chrom_size, args.windowsize + pos])}\t{avg_idat_afr}\t{avg_idat_eur}\t{avg_idat_eas}\n')
            pos += args.stride
    outfile.close()


def get_standardized_idat_scores_regions_helper(args):
    regions, df, local_ancestry_dict, dist_step_size, min_cov_dat = args
    iDAT_AFR = []
    iDAT_EUR = []
    iDAT_EAS = []
    for chrom, start, end in zip(regions.chrom.values, regions.start.values, regions.end.values):
        local_ancestry = local_ancestry_dict[chrom]
        # TODO: SHOULD I CALCULATE IT FOR EVERY POSITION IN A GIVEN REGION RATHER THAN JUST THE MIDPOINT?
        _, idat_scores = calculate_idat_score((local_ancestry, [int((end + start) / 2)], dist_step_size, min_cov_dat))
        iDAT_AFR.append(idat_scores[0][0])
        iDAT_EUR.append(idat_scores[0][1])
        iDAT_EAS.append(idat_scores[0][2])
    # standardize to 0 mean and 1 standard deviation
    regions['standardized_iDAT_AFR'] = (np.array(iDAT_AFR) - df.idat_AFR.mean()) / df.idat_AFR.std()
    regions['standardized_iDAT_EUR'] = (np.array(iDAT_EUR) - df.idat_EUR.mean()) / df.idat_EUR.std()
    regions['standardized_iDAT_EAS'] = (np.array(iDAT_EAS) - df.idat_EAS.mean()) / df.idat_EAS.std()
    return regions


def get_standardized_idat_scores_regions(args):
    df = pd.concat([pd.read_csv(idat, sep='\t', names=['chrom', 'pos', 'idat_AFR', 'idat_EUR', 'idat_EAS'])
                    for idat in args.idat_scores])
    regions = pd.read_csv(args.regions, sep='\t', header=0)
    if not 'chrom' in regions.columns.values:
        regions = pd.read_csv(args.regions, sep='\t', usecols=[0, 1, 2], names=['chrom', 'start', 'end'])

    local_ancestry_pattern = re.sub("chr[0-9]+", "chr{chrom}", args.local_ancestry_pattern)
    local_ancestry_pattern = re.sub("phase[0, 1]", "phase{phase}", local_ancestry_pattern)
    local_ancestry_dict = {}
    for chrom in regions.chrom.unique():
        local_ancestry_dict[chrom] = pd.concat([pd.read_csv(local_ancestry_pattern.format(chrom=chrom.replace('chr', ''), phase='0'),
                                                            sep='\t', header=0),
                                                pd.read_csv(local_ancestry_pattern.format(chrom=chrom.replace('chr', ''), phase='1'),
                                                            sep='\t', header=0)])
    splits = np.array_split(regions, regions.shape[0] // 8 + 1)
    # calculate iDAT for each region
    ready_to_map = [(split, df, local_ancestry_dict, args.dist_step_size, args.min_cov_dat)
                    for split in splits]
    pool = mp.Pool(processes=args.threads)
    results = pool.map(get_standardized_idat_scores_regions_helper, ready_to_map)
    pool.close()
    pool.join()
    regions = pd.concat(results).sort_values(['chrom', 'start'])
    regions.to_csv(args.output, sep='\t', header=True, index=False)


def main(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    subparser = subparsers.add_parser('idat_scores')
    subparser.add_argument('-l', '--local_ancestry', nargs='+', help='Local ancestry tracks in both phases inferred '
                                                                     'with FLARE. Must contain tracks labeled as '
                                                                     'AFR, EUR, and EAS')
    subparser.add_argument('-g', '--genomefile', help='Path to genomefile as expected by BEDTools')
    subparser.add_argument('-o', '--output', help='Path to output BED file containg iDAT scores for all '
                                                  'positions.')
    subparser.add_argument('-n', '--number_total_sites', type=int, help='Number of total random genomic sites for which'
                                                                        ' to calculate standardized iDAT scores. '
                                                                        '[10000]', default=10000)
    subparser.add_argument('-c', '--chrom', type=int, help='Chromosome')
    subparser.add_argument('-m', '--min_cov_dat', type=float, help='Minimal DAT score considered when calculating AUC. '
                                                                   '[0.1]', default=0.1)
    subparser.add_argument('-d', '--dist_step_size', type=int, help='Step size with which to move away from reference '
                                                                    'point when calculating DAT scores. [10000]',
                           default=10000)
    subparser.add_argument('-t', '--threads', type=int, help='Number of threads')
    subparser.set_defaults(func=get_idat_scores)

    subparser = subparsers.add_parser('standardized_idat_scores')
    subparser.add_argument('-i', '--idat_scores', nargs='+', help='Paths to unstandardized genome-wide iDAT scores.')
    subparser.add_argument('-g', '--genomefile', help='Path to genomefile as expected by BEDTools')
    subparser.add_argument('-o', '--output', help='Path to output BED file containg standardized iDAT scores for all '
                                                  'windows.')
    subparser.add_argument('-w', '--windowsize', type=int, help='Window size for which to calculate average '
                                                                'standardized iDAT score. [20 Mb]', default=20e6)
    subparser.add_argument('-s', '--stride', type=int, help='Steps size with which to move window. [1 Mb]', default=1e6)
    subparser.set_defaults(func=get_standardized_idat_scores_windows)

    subparser = subparsers.add_parser('regions')
    subparser.add_argument('-r', '--regions', help='BED file with header specifying regions for which to calculate '
                                                   'standardized iDAT scores')
    subparser.add_argument('-l', '--local_ancestry_pattern', help='Example file name pattern of local ancestry tracks in both '
                                                                  'phases inferred with FLARE. Must contain tracks '
                                                                  'labeled as AFR, EUR, and EAS. Must contain '
                                                                  'substrings matching "chr[0-9]+" and "phase[0, 1]".')
    subparser.add_argument('-i', '--idat_scores', nargs='+', help='Paths to genome-wide unstandardized iDAT scores.')
    subparser.add_argument('-d', '--dist_step_size', type=int, help='Step size with which to move away from reference '
                                                                    'point when calculating DAT scores. [1000]',
                           default=1000)
    subparser.add_argument('-m', '--min_cov_dat', type=float, help='Minimal DAT score considered when calculating AUC. '
                                                                   '[0.1]', default=0.1)
    subparser.add_argument('-o', '--output', help='Path to output BED file containing standardized iDAT scores for '
                                                  'each region.')
    subparser.add_argument('-t', '--threads', type=int, help='Number of threads. [1]', default=1)
    subparser.set_defaults(func=get_standardized_idat_scores_regions)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
