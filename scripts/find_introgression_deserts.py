#!/usr/bin/env python3
import pandas as pd
import numpy as np
from pybedtools import BedTool, helpers
import multiprocessing as mp
import sys
import argparse
from identify_putatively_selected_introgressed_segments_in_african_americans import round_local_ancestry_haplotypes, \
    calculate_expectation_and_probability_helper, sample_control_sequences_helper, get_windows, get_overlap_mask
import logging
from multiprocessing_logging import install_mp_handler


def get_introgressed_bp_per_window(args):
    """
    Annotate windows with number bp that are masked and introgressed bp across all individuals
    :param args: tuple, (introgression_coverage, genomefile, window_size, step size, masked, max_masked)
    :return: (pd.DataFrame, pd.DataFrame), genomic windows with information about fraction of masked bp and
    introgressed bp across all individuals as well as a dataframe of remaining callable windows
    """
    introgression_coverage, genomefile, window_size, stride, masked, max_masked, n_indvs = args
    # get windows
    windows = get_windows(genomefile, window_size, stride)
    # aggregate masked bp per window
    windows = BedTool.from_dataframe(windows).intersect(BedTool.from_dataframe(masked),
                                                        wao=True).groupby(g=[1, 2, 3],
                                                                          c=7, o=['sum']).to_dataframe(names=['chrom',
                                                                                                              'start',
                                                                                                              'end',
                                                                                                              'masked'])
    # filter out windows in which more than half of the bases are masked
    callable_windows = windows.loc[windows.masked / window_size < (1 - max_masked)]
    if callable_windows.shape[0] == 0:
        return (None, None)
    else:
        introgressed = BedTool.from_dataframe(callable_windows).intersect(BedTool.from_dataframe(introgression_coverage),
                                                                          wao=True).groupby(g=[1, 2, 3], c=9,
                                                                                            o=['sum']).to_dataframe(names=['chrom',
                                                                                                                           'start',
                                                                                                                           'end',
                                                                                                                           'introgressed'])

        callable_windows['introgressed'] = np.where(introgressed.introgressed.values < 0, 0, introgressed.introgressed.values)
        callable_windows.introgressed /= n_indvs
        callable_windows.introgressed /= window_size
        # define significant windows
        introgression_cutoff = np.percentile(callable_windows.introgressed.values, 1)
        deserts = BedTool.from_dataframe(callable_windows[callable_windows.introgressed <= introgression_cutoff]).merge().to_dataframe()
        return (deserts, callable_windows)


def get_control_sequences_for_deserts(callable_windows, deserts, output, genomefile, reps, threads):
    # get random control sequences
    if deserts.shape[0] > 0:
        callable_windows = BedTool.from_dataframe(callable_windows).subtract(BedTool.from_dataframe(deserts)).merge().to_dataframe(names=['chrom', 'start', 'end'])
        callable_windows['length'] = callable_windows.end - callable_windows.start
        deserts['length'] = deserts.end - deserts.start
        ready_to_map_random = [(deserts, callable_windows, genomefile, np.random.randint(0, 1e7)) for n in range(reps)]
        pool = mp.Pool(processes=threads)
        random_control_segments_replicates = pool.map(sample_control_sequences_helper, ready_to_map_random)
        pool.close()
        pool.join()
        for i, random_control_segments in enumerate(random_control_segments_replicates):
            # random_control_segments.chrom = random_control_segments.chrom.astype(int)
            random_control_segments.start = random_control_segments.start.astype(int)
            random_control_segments.end = random_control_segments.end.astype(int)
            random_control_segments.to_csv(f'{output}_{i}.bed', index=False, header=True, sep='\t')
    else:
        for i in range(reps):
            random_control_segments = pd.DataFrame(columns=['chrom', 'start', 'end'])
            random_control_segments.to_csv(f'{output}_{i}.bed', index=False, header=True, sep='\t')


def add_windowsize_specific_rank_point(callable_windows):
    """
    At window sizes specific rank of a window in terms of introgression frequency
    :param callable_windows: pd.DataFrame, callable windows with inferred introgression frequencies
    :return: pd.DataFrame, callable_windows with rank
    """
    # infer window sizes
    callable_windows['windowsize'] = np.round((callable_windows.end - callable_windows.start) / 1e6)
    callable_windows.windowsize = callable_windows.windowsize.astype(int)
    # infer rank
    for ws in callable_windows.windowsize.unique():
        callable_windows.loc[callable_windows.windowsize == ws,
                             'rank'] = callable_windows.loc[callable_windows.windowsize == ws,
                                                            'introgressed'].rank(pct=True)
    return callable_windows


def get_local_ancestry_for_each_region_helper(args):
    merged_callable_windows, local_anc_fn_pattern, chrom, columns = args
    # get callable windows on current chromosome
    c_callable = merged_callable_windows.loc[merged_callable_windows.chrom == chrom,
                                             ['chrom', 'start', 'end']]
    # intersect with local ancestry file phase0
    if not isinstance(chrom, str):
        chrom = f"chr{chrom}"
    local_anc_freq0 = BedTool.from_dataframe(c_callable).intersect(
        local_anc_fn_pattern.format(chrom=chrom.replace('chr', ''), phase=0),
        wao=True).to_dataframe(names=columns)
    local_anc_freq0.drop(['chrom_la', 'start_la', 'end_la', 'IID_la'], axis=1, inplace=True)
    # group local ancestry by region
    local_anc_freq0 = local_anc_freq0.groupby(['chrom', 'start', 'end', 'la'])['overlap'].sum()

    # intersect with local ancestry file phase1
    local_anc_freq1 = BedTool.from_dataframe(c_callable).intersect(
        local_anc_fn_pattern.format(chrom=chrom.replace('chr', ''), phase=1),
        wao=True).to_dataframe(names=columns)
    local_anc_freq1.drop(['chrom_la', 'start_la', 'end_la', 'IID_la'], axis=1, inplace=True)
    # group local ancestry by region
    local_anc_freq1 = local_anc_freq1.groupby(['chrom', 'start', 'end', 'la'])['overlap'].sum()

    # merge both phases
    local_anc_freq_merged = pd.DataFrame(local_anc_freq0).join(local_anc_freq1, rsuffix='_1')
    local_anc_freq_merged.overlap += local_anc_freq_merged.overlap_1
    local_anc_freq_merged.drop('overlap_1', axis=1, inplace=True)
    local_anc_freq_merged.reset_index(inplace=True)

    # add local ancestry info to callable windows
    local_anc_freq_merged.set_index(['chrom', 'start', 'end'], inplace=True)
    c_callable.set_index(['chrom', 'start', 'end'], inplace=True)
    for label in local_anc_freq_merged.la.unique():
        c_callable = c_callable.join(
            local_anc_freq_merged[local_anc_freq_merged.la == label])
        c_callable.drop('la', axis=1, inplace=True)
        c_callable.rename({'overlap': label}, axis=1, inplace=True)
    # calculate local ancestry proportions
    anc_sum = c_callable.sum(axis=1)
    for label in local_anc_freq_merged.la.unique():
        c_callable[label] /= anc_sum
    return c_callable


def get_local_ancestry_for_each_region(merged_callable_windows, local_anc_fn_pattern, threads):
    columns = ['chrom', 'start', 'end']
    columns.extend(['chrom_la', 'start_la', 'end_la', 'la', 'IID_la'])
    columns.append('overlap')
    ready_to_map = [(merged_callable_windows, local_anc_fn_pattern, chrom, columns)
                    for chrom in merged_callable_windows.chrom.unique()]
    # get_local_ancestry_for_each_region_helper(ready_to_map[0])
    pool = mp.Pool(processes=threads)
    local_anc_per_chrom = pool.map(get_local_ancestry_for_each_region_helper, ready_to_map)
    pool.close()
    pool.join()
    # concat all chromosomes
    local_anc_per_chrom = pd.concat(local_anc_per_chrom).sort_values(['chrom', 'start', 'end'])
    merged_callable_windows = merged_callable_windows.set_index(['chrom', 'start', 'end']).join(
        local_anc_per_chrom).reset_index()
    return merged_callable_windows


def identify_new_deserts_in_admixed_population(callable_windows, local_anc_fn_pattern, callable_windows_references,
                                               callable_windows_references_labels, n_indvs,  n_eur_indv,
                                               n_eas_indv, threads):
    # merge callable windows and introgression frequency with those from references
    merged_callable_windows = callable_windows[callable_windows['rank']
                                               < 0.05].set_index(['chrom', 'start', 'end']).join(
    callable_windows_references[0].set_index(['chrom', 'start', 'end']),
        rsuffix=f'_{callable_windows_references_labels[0]}')
    for cw, label in zip(callable_windows_references[1:], callable_windows_references_labels[1:]):
        merged_callable_windows = merged_callable_windows.join(cw.set_index(['chrom', 'start', 'end']),
                                                               rsuffix=f'_{label}')
    merged_callable_windows.reset_index(inplace=True)
    # add local ancestry information
    merged_callable_windows = get_local_ancestry_for_each_region(merged_callable_windows, local_anc_fn_pattern, threads)
    # format column names to use existing model
    merged_callable_windows.rename({'introgressed_EUR': 'freq_eur',
                                    'introgressed_EAS': 'freq_eas', 'introgressed': 'freq',
                                    'AFR': 'freq_la_afr', 'EUR': 'freq_la_eur',
                                    'EAS': 'freq_la_eas'}, axis=1, inplace=True)

    # get number of local ancestry haplotypes pre region
    merged_callable_windows.freq_la_afr = np.floor(merged_callable_windows.freq_la_afr * n_indvs)
    merged_callable_windows.freq_la_eur = np.floor(merged_callable_windows.freq_la_eur * n_indvs)
    merged_callable_windows.freq_la_eas = np.floor(merged_callable_windows.freq_la_eas * n_indvs)
    merged_callable_windows = round_local_ancestry_haplotypes(merged_callable_windows, n_indvs)
    # get number of introgressed haplotypes per region
    merged_callable_windows.freq = np.ceil(merged_callable_windows.freq * n_indvs)
    # avoid negatives --> that is there was 0 introgression (bedtools returns -1)
    merged_callable_windows.freq_eur = np.where(merged_callable_windows.freq_eur.values < 0, 0,
                                                merged_callable_windows.freq_eur.values)
    merged_callable_windows.freq_eas = np.where(merged_callable_windows.freq_eas.values < 0, 0,
                                                merged_callable_windows.freq_eas.values)
    merged_callable_windows.freq = np.where(merged_callable_windows.freq.values < 0, 0,
                                            merged_callable_windows.freq.values)
    splits = np.array_split(merged_callable_windows.sample(frac=1,  random_state=42), merged_callable_windows.shape[0] // 8)
    pool = mp.Pool(processes=threads)
    ready_to_map = [(split, n_eur_indv, n_eas_indv, n_indvs, True, i, len(splits)) for i, split in enumerate(splits)]
    dfs = pool.map(calculate_expectation_and_probability_helper, ready_to_map)
    logging.info('Finished all splits')
    pool.close()
    pool.join()
    df = pd.concat(dfs)
    df.sort_values(['chrom', 'start'], inplace=True)
    n_independent_regions = BedTool.from_dataframe(df).merge().to_dataframe().shape[0]
    df['corrected'] = np.where(df.probabilities.values * n_independent_regions > 1.0, 1.0,
                               df.probabilities.values * n_independent_regions)
    df_pvalues = df.copy()
    # apply ancestry filter and require some introgression in Europe
    df = df[(df.freq < df.expectations) &
            (df.corrected < 0.05) & (df.freq_eur > 0) &
            (df.freq_la_afr / n_indvs >= 0.5) &
            (df.freq_la_eur / n_indvs >= 0.1) &
            (df.freq_la_eas / n_indvs < 0.05)]
    df = df[(df.freq < df.expectations) & (df.corrected < 0.05)]
    df = BedTool.from_dataframe(df).merge().to_dataframe()
    return df, df_pvalues


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--introgression_coverage', help='File with introgression coverage in bp across all '
                                                         'individuals for each genomic window')
    parser.add_argument('--introgression_coverage_masked', help='File with introgression coverage in bp across all '
                                                                'individuals for each genomic window '
                                                                'after masking segments also found in Africa')
    parser.add_argument('--ibdmix', help='IBDmix calls for all individuals')
    parser.add_argument('--ibdmix_masked', help='IBDmix calls for all individuals after masking segments '
                                                'also found in Africa')
    parser.add_argument('--mask_pattern', help='String pattern of mask files, must allow string formatting of "chrom"')
    parser.add_argument('--genomefile', help='Path to genomefile as BEDTools expects it')
    parser.add_argument('-o', '--output', help='Output file with introgression deserts')
    parser.add_argument('-or', '--output_random', help='Prefix for output files that will contain random '
                                                       'control segments. "_{n}.bed" suffix will be appended',
                        default=None)
    parser.add_argument('-orn', '--output_random_new', help='Prefix for output files that will contain random '
                                                            'control segments for new desert like regions. '
                                                            '"_{n}.bed" suffix will be appended',
                        default=None)
    parser.add_argument('--output_ranks', help='File with number of masked and introgressed bases per window '
                                               'as well as the windowsize specific rank associated with number '
                                               'of introgressed bases.')
    parser.add_argument('--output_ranks_masked', help='File with number of masked and introgressed bases per window '
                                                      'as well as the windowsize specific rank associated with number '
                                                      'of introgressed bases after masking introgressed '
                                                      'segments also found in Africa.')
    parser.add_argument('--reps', type=int, default=100,
                        help='Number of bootstrap sets to create from control segments. [100]')
    parser.add_argument('--window_sizes', nargs='+', type=int,
                        help='Sizes of windows in bp genome is to be split into '
                             '[8e6, 9e6, 10e6, 11e6, 12e6, 13e6, 14e6, 15e6]',
                        default=[8e6, 9e6, 10e6, 11e6, 12e6, 13e6, 14e6, 15e6])
    parser.add_argument('--stride', type=int, help='Step size of windows in bp, [100000]', default=100000)
    parser.add_argument('--max_masked', type=float, help='Maximal fraction of bp that are allowed to be masked in '
                                                         'a window to be considered [0.5]', default=0.5)
    parser.add_argument('--chrom', help='Consider only a specific chromosome. [None]', default=None)
    parser.add_argument('--super_pop', help='Population to consider')
    parser.add_argument('--threads', type=int, help='Number of CPUS [8]', default=8)
    parser.add_argument('-l', '--lai', help='Local ancestry filename pattern. Must allow string formatting with '
                                            'chrom and phase. Only needed when looking for novel '
                                            'deserts in admixed population. [None]', default=None)
    parser.add_argument('--callable_windows_references', help='Paths to dataframes with introgression frequencies for '
                                                              'callable windows for different reference populations. '
                                                              'Only needed when looking for novel deserts in '
                                                              'admixed population. [None]', default=None, nargs='+')
    parser.add_argument('--callable_windows_references_labels', help='Population labels for above specified files '
                                                                     'Only needed when looking for novel deserts in '
                                                                     'admixed population. [None]', default=None,
                        nargs='+')
    parser.add_argument('--unique_deserts', help='Filename for file with unique deserts.')
    parser.add_argument('--unique_deserts_pvalues', help='Filename for file with unique deserts with pvalues.')
    parser.add_argument('--tmp_dir', help='Temporary directory for pybedtools. Default=/tmp/', default='/tmp/')

    args = parser.parse_args()
    logging.basicConfig(filename='finding_introgression_deserts.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')

    logger = logging.getLogger()
    install_mp_handler()
    if args.chrom:
        mask = pd.read_csv(args.mask_pattern.format(chrom=args.chrom), sep='\t', names=['chrom', 'start', 'end'])

    else:
        try:
            mask = pd.concat([pd.read_csv(args.mask_pattern.format(chrom=c_chrom), sep='\t',
                                          names=['chrom', 'start', 'end'])
                              for c_chrom in range(1, 23)])
        # simulations only have 10 chromosomes
        except FileNotFoundError:
            mask = pd.concat([pd.read_csv(args.mask_pattern.format(chrom=c_chrom), sep='\t',
                                          names=['chrom', 'start', 'end'])
                              for c_chrom in range(1, 11)])
    helpers.set_tempdir(args.tmp_dir)
    mask.chrom = [f'chr{chrom}' for chrom in mask.chrom]
    mask.sort_values(['chrom', 'start', 'end'], inplace=True)
    mask.reset_index(inplace=True, drop=True)
    # get masked base pair for each non overlapping window in parallel --> can aggregate to window sizes later
    windows = get_windows(args.genomefile, args.stride, args.stride)
    ready_to_map_mask = [(windows, mask, chrom) for chrom in windows.chrom.unique()]
    pool = mp.Pool(processes=args.threads)
    masked = pool.map(get_overlap_mask, ready_to_map_mask)
    pool.close()
    pool.join()
    masked = pd.concat(masked)
    masked.sort_values(['chrom', 'start'])
    """
    Identify traditional introgression deserts
    """
    # load introgression_coverage
    introgression_coverage = pd.read_csv(args.introgression_coverage, sep='\t',
                                         names=['chrom', 'start', 'end', 'super_pop', 'coverage'])
    # load ibdmix
    ibdmix = pd.read_csv(args.ibdmix, sep='\t', names=['chrom', 'start', 'end', "LOD", "IID", "pop", "super_pop"])
    n_indvs = ibdmix.loc[ibdmix.super_pop == args.super_pop, 'IID'].unique().shape[0]
    ready_to_map_introgression = [(introgression_coverage, args.genomefile, w, args.stride, masked,
                                   args.max_masked, n_indvs) for w in args.window_sizes]
    # get_introgressed_bp_per_window(ready_to_map_introgression[0])
    # get_masked_and_introgressed_bp_per_window(ready_to_map[0])
    pool = mp.Pool(processes=args.threads)
    results = pool.map(get_introgressed_bp_per_window, ready_to_map_introgression)
    pool.close()
    pool.join()

    deserts = [res[0] for res in results]
    deserts = pd.concat(deserts)
    # merge deserts across all window sizes
    if deserts.shape[0] > 0:
        deserts.sort_values(['chrom', 'start'], inplace=True)
        deserts = BedTool.from_dataframe(deserts).merge().to_dataframe(names=['chrom', 'start', 'end'])
    deserts.to_csv(args.output, sep='\t', header=False, index=False)
    callable_windows = pd.concat([res[1] for res in results])
    callable_windows.sort_values(['chrom', 'start'], inplace=True)
    callable_windows = add_windowsize_specific_rank_point(callable_windows)
    callable_windows.to_csv(args.output_ranks, sep='\t', header=True, index=False)
    if args.output_random:
        # sample control segments
        # get windows that are not in the lowest 99th percentile
        callable_windows = callable_windows[callable_windows['rank'] > 0.01]
        callable_windows = BedTool.from_dataframe(callable_windows).merge().to_dataframe(names=['chrom', 'start', 'end'])
        get_control_sequences_for_deserts(callable_windows, deserts, args.output_random, args.genomefile, args.reps,
                                          args.threads)
    """
    Identify novel introgression using IBDmix calls after applying the "African mask"
    """
    # load introgression_coverage
    introgression_coverage_masked = pd.read_csv(args.introgression_coverage_masked, sep='\t',
                                                names=['chrom', 'start', 'end', 'super_pop', 'coverage'])
    # load ibdmix
    ibdmix_masked = pd.read_csv(args.ibdmix_masked, sep='\t', names=['chrom', 'start', 'end', "LOD", "IID", "pop", "super_pop"])
    n_indvs = ibdmix_masked.loc[ibdmix_masked.super_pop == args.super_pop, 'IID'].unique().shape[0]
    ready_to_map_introgression = [(introgression_coverage_masked, args.genomefile, w, args.stride, masked,
                                   args.max_masked, n_indvs) for w in args.window_sizes]
    # get_masked_and_introgressed_bp_per_window(ready_to_map[0])
    pool = mp.Pool(processes=args.threads)
    results = pool.map(get_introgressed_bp_per_window, ready_to_map_introgression)
    pool.close()
    pool.join()

    callable_windows_masked = pd.concat([res[1] for res in results])
    callable_windows_masked.sort_values(['chrom', 'start'], inplace=True)
    callable_windows_masked = add_windowsize_specific_rank_point(callable_windows_masked)
    callable_windows_masked.to_csv(args.output_ranks_masked, sep='\t', header=True, index=False)
    # identify novel desert-like regions in admixed population
    if not args.lai is None:
        assert args.callable_windows_references is not None and args.callable_windows_references_labels is not None, \
            "If local ancestry information is provided, callable_windows_references and " \
            "callable_windows_reference_labels must also be specified"
        callable_windows_references = [pd.read_csv(ref, sep='\t', header=0) for ref in args.callable_windows_references]
        callable_windows_references_labels = args.callable_windows_references_labels
        n_eur_indv = ibdmix.loc[ibdmix.super_pop == 'EUR', 'IID'].unique().shape[0]
        n_eas_indv = ibdmix.loc[ibdmix.super_pop == 'EAS', 'IID'].unique().shape[0]
        new_deserts, new_deserts_pvalues = identify_new_deserts_in_admixed_population(callable_windows_masked.copy(),
                                                                                      args.lai,
                                                                                      callable_windows_references,
                                                                                      callable_windows_references_labels,
                                                                                      n_indvs, n_eur_indv,
                                                                                      n_eas_indv, args.threads)
        new_deserts.to_csv(args.unique_deserts, sep='\t', header=False, index=False)
        new_deserts_pvalues.to_csv(args.unique_deserts_pvalues, sep='\t', header=True, index=False)
        if args.output_random_new:
            # sample control segments
            # get windows that are not in the lowest 99th percentile
            callable_windows_potential_new_deserts = callable_windows_masked[callable_windows_masked['rank'] > 0.05].copy()
            callable_windows_potential_new_deserts = BedTool.from_dataframe(
                callable_windows_potential_new_deserts).merge().to_dataframe(names=['chrom', 'start', 'end'])
            get_control_sequences_for_deserts(callable_windows_potential_new_deserts, new_deserts,
                                              args.output_random_new, args.genomefile, args.reps, args.threads)


if __name__ == '__main__':
    main(sys.argv[1:])
