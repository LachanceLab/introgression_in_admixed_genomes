#!/usr/bin/env python
import sys
import argparse
import pandas as pd
import numpy as np
from scipy.stats import binom, norm, hmean, combine_pvalues, t
from pybedtools import BedTool, helpers
from sklearn.metrics import pairwise_distances
import multiprocessing as mp
import logging
from multiprocessing_logging import install_mp_handler


def get_windows(genomefile, window_size, stride=100000):
    """
    Split genome into windows of defined size with defined step size
    :param genomefile: str, path to genomefile as BEDTools would expect it
    :param window_size: int, window size
    :param stride: int, step size
    :return: pd.DataFrame, genomic windows
    """
    windows = BedTool().window_maker(g=genomefile, s=stride,
                                     w=window_size).to_dataframe()
    # drop short windows towards the end of chromosome that are already covered
    chrom_end_windows = windows[windows.end - windows.start < window_size]
    windows.drop(chrom_end_windows[chrom_end_windows.duplicated('chrom')].index, inplace=True)
    return windows


def get_overlap_mask(args):
    windows, mask, chrom = args
    c_windows = windows[windows.chrom == chrom]
    c_mask = mask[mask.chrom == chrom]
    if c_mask.shape[0] == 0:
        mask_tmp = mask.copy()
        mask_tmp.chrom = [int(chrom.replace('chr', '')) for chrom in mask_tmp.chrom]
        c_mask = mask_tmp[mask_tmp.chrom == chrom]

    masked = BedTool.from_dataframe(c_windows).intersect(BedTool.from_dataframe(c_mask),
                                                         wao=True).groupby(g=[1, 2, 3],
                                                                           c=7, o=['sum']).to_dataframe(names=['chrom',
                                                                                                               'start',
                                                                                                               'end',
                                                                                                               'masked'])
    return masked


def meff_gao(x):
    """
    Estimate number of effective tests using method by Gao et al. (https://doi.org/10.1002/gepi.20310)
    :param x: array-like, eigenvalues of correlation matrix
    :return: int, number of effective tests
    """
    return np.where(np.cumsum(np.arange(1, x.shape[0] + 1) * x) /
                    sum(np.arange(1, x.shape[0] + 1) * x) < 0.995)[0].shape[0]


def meff_li_ji(x):
    """
    Estimate number of effective tests using method by Li & Ji (https://doi.org/10.1038/sj.hdy.6800717)
    :param x: array-like, eigenvalues of correlation matrix
    :return: int, number of effective tests
    """
    return np.abs(np.where(x >= 1, 1, 0) + (x - np.floor(x))).sum()


def meff_nyholt(x):
    """
    Estimate number of effective tests using method by Nyholt (https://doi.org/10.1086/383251)
    :param x: array-like, eigenvalues of correlation matrix
    :return: int, number of effective tests
    """
    return 1 + (x.shape[0] - 1) * (1 - (np.std(x) ** 2) / x.shape[0])


def meff_galeway(x):
    """
    Estimate number of effective tests using method by Galeway (https://doi.org/10.1002/gepi.20408)
    :param x: array-like, eigenvalues of correlation matrix
    :return: int, number of effective tests
    """
    return np.floor(np.sum(np.sqrt(np.where(x < 0, 0, x))) ** 2 / np.sum(x))


def parse_genome_file(file):
    """
    Parse BedTools genome file to get chromosome sizes
    :param file: str, path to genome file for BedTools
    :return: dict, chromosome sizes
    """
    chrom_sizes = {}
    with open(file, 'r') as genome:
        for line in genome:
            line = line.strip().split('\t')
            chrom_sizes[line[0]] = int(line[1])
    genome.close()
    return chrom_sizes


def ecdf(x):
    '''
    Modified version from statsmodels.stats.multitest
    (https://www.statsmodels.org/0.8.0/_modules/statsmodels/stats/multitest.html)

    no frills empirical cdf used in fdrcorrection
    '''
    nobs = len(x)
    return (np.arange(1, nobs + 1) - 1) / (float(nobs) - 1)


def fdrcorrection_meff(pvals, meff, alpha=0.05, is_sorted=False):
    """
    Modified version from statsmodels.stats.multitest
    (https://www.statsmodels.org/0.8.0/_modules/statsmodels/stats/multitest.html)

    FDR-procedure according to Benjamini-Hochberg but using the effective number of independent tests as proposed by
    Li and Ji (https://www.nature.com/articles/6800717)
    :param pvals: array-like, set of p-values of the individual tests
    :param meff: int, estimated effective number of independent tests
    :param alpha: float, error rate
    :param is_sorted: boolean, whether pvalues are sorted or not
    :return: array, array, rejected (boolean array indicating if hypothesis is rejected or not),
                           pvalues (pvalues adjusted for multiple testing to limit FDR)
    """
    pvals = np.asarray(pvals)
    # sort pvalues
    if not is_sorted:
        pvals_sortind = np.argsort(pvals)
        pvals_sorted = np.take(pvals, pvals_sortind)
    # already sorted
    else:
        pvals_sorted = pvals  # alias

    # calculate empirical cdf according to Li and Ji
    ecdffactor = ecdf(pvals_sorted)

    # perform hypothesis testing according to Li and Ji
    reject = pvals_sorted < alpha / meff + ecdffactor * (alpha - alpha / meff)
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
        reject[:rejectmax] = True
    # adjust p-values for multiple testing
    # solve p <= q/meff + ecdffactor * (q - q/meff) for q
    pvals_corrected_raw = pvals_sorted * meff / (ecdffactor * (meff - 1) + 1)
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    del pvals_corrected_raw
    pvals_corrected[pvals_corrected > 1] = 1
    if not is_sorted:
        pvals_corrected_ = np.empty_like(pvals_corrected)
        pvals_corrected_[pvals_sortind] = pvals_corrected
        del pvals_corrected
        reject_ = np.empty_like(reject)
        reject_[pvals_sortind] = reject
        return reject_, pvals_corrected_
    else:
        return reject, pvals_corrected


def add_count(x, c):
    x[int(x[-1])] += c
    return x


def sample_cols(x):
    return np.random.choice([0, 1, 2], size=1, p=x)


def estimate_selection_coefficient(p0, p1, t=15):
    return (np.log((p1 * (p0 - 1)) / (p0 * (p1 - 1)))) / t


def neanderthal_frequency_dist_with_uncertainty(p, n, nstd=3.29, discretize=False):
    """
    Compute Neanderthal introgression frequency distribution for reference population
    :param p: float, measured allele frequency
    :param n: int, sample size (individuals)
    :param nstd: int, number of standard deviation out to consider
    :return: initialized normal distribution, lower bound (float), upper bound (float)
    """
    # calculate p' according to agresti-coull --> don't do the correction when p == 0
    # if p > 0:
    p_adj = (p * n + 2) / (n + 4)
    # else:
    #     p_adj = p
    # define std of allele frequency dist according to agresti-coull
    std = np.sqrt((p_adj * (1 - p_adj)) / (n + 4))
    # initialize pdf
    pdf = norm(loc=p_adj, scale=std)
    # get CI
    lower_bound = p_adj - nstd * std
    upper_bound = p_adj + nstd * std
    possible_vals = np.arange(max([lower_bound, 0]), min([1, upper_bound]) + 0.001, 0.001)
    if possible_vals[-1] < upper_bound:
        possible_vals = np.concatenate([possible_vals, [upper_bound]])
    pdf_vals = pdf.pdf(possible_vals)
    # define lower bound --> cannot be lower than 0
    if lower_bound < 0.0:
        pdf_vals[0] += pdf.pdf(np.arange(lower_bound, 0, 0.001)).sum()
        lower_bound = 0.0
    # define upper bound
    if upper_bound > 1.0:
        upper_bound = 1.0
        pdf_vals[-1] += pdf.pdf(np.arange(1, upper_bound, 0.001)).sum()
    # discretize
    if discretize:
        discrete_lower = np.floor(lower_bound * n)
        discrete_upper = np.ceil(upper_bound * n)
        possible_vals = np.arange(discrete_lower, discrete_upper + 1, 1)
        pdf_vals_discrete = pdf.pdf(possible_vals / n)
        # collapse out of bounds
        if lower_bound == 0:
            pdf_vals_discrete[0] = pdf_vals[0]
        if upper_bound == 1.0:
            pdf_vals_discrete[-1] = pdf_vals[-1]
        pdf_vals = pdf_vals_discrete / pdf_vals_discrete.sum()
    else:
        pdf_vals /= pdf_vals.sum()

    return pdf, lower_bound, upper_bound, possible_vals, pdf_vals


def get_min_and_max_contributions(trials, p_lower, p_upper, tol=1e-20):
    """
    Only consider contribution that have a probability >= tolerance to limit search space
    :param trials: int, number of bernoulli trials
    :param p_lower: float, lower success rate
    :param p_upper: float, upper success rate
    :param tol: float, min probability in order to include a value in the range
    :return int, int; min and max contributions
    """
    prob_contr = binom.pmf(np.arange(0, trials + 1), trials, p_lower)
    min_contr = np.where(prob_contr >= tol)[0][0]
    prob_contr = binom.pmf(np.arange(0, trials + 1), trials, p_upper)
    max_contr = np.where(prob_contr >= tol)[0][-1] + 1
    return min_contr, max_contr


# def get_dimensions(afr, eur, eas, p_amr, comparison):
def get_dimensions(eur, eas, p_amr, comparison):
    """
    Get indices of relevant combinations of contributions from each ancestry
    :param eur: np.ndarray, possible EUR contributions
    :param eas: np.ndarray, possible EAS contributions
    :param p_amr: np.ndarray, AMR observation
    :param comparison: str, 'smaller' or 'greater', indicating whether observation in AMR is larger
                       or smaller than expected
    :return: tuple, (array, array, array) with relevant indices across all 3 axes
    """
    # depletion
    if comparison == 'smaller':
        # dimensions = np.where(afr + eur + eas <= p_amr)
        dimensions = np.where(eur + eas <= p_amr)
    # enrichment
    else:
        # dimensions = np.where(afr + eur + eas >= p_amr)
        dimensions = np.where(eur + eas >= p_amr)
    return dimensions


def calculate_expectation_and_probability_helper(args):
    """
    Calculate expectation and probability of observed introgressed frequency in target population given local ancestry
    frequencies and introgression frequencies in reference populations
    """
    # df, n_eur_indv, n_afr_indv, n_eas_indv, n_amr_indv, calculate_probabilities, c_split, n_splits = args
    df, n_eur_indv, n_eas_indv, n_amr_indv, calculate_probabilities, c_split, n_splits = args

    logging.info(f"Working on split {c_split}/{n_splits}\n")
    expectations = []
    probabilities = []
    selection_coefficients = []
    corrected_n_haplotypes = []
    # selection_coefficients_low_ci = []
    # selection_coefficients_high_ci = []
    # for n_haplotypes, trials_afr, trials_eur, trials_eas, p_afr, p_eur, p_eas in df.loc[:, ['freq', 'freq_la_afr',
    #                                                                                         'freq_la_eur',
    #                                                                                         'freq_la_eas',
    #                                                                                         'freq_afr', 'freq_eur',
    #                                                                                         'freq_eas']].values:
    for n_haplotypes, trails_afr, trials_eur, trials_eas,  p_eur, p_eas in df.loc[:, ['freq', 'freq_la_afr',
                                                                                      'freq_la_eur', 'freq_la_eas',
                                                                                      'freq_eur', 'freq_eas']].values:

        p_eur_pdf, p_eur_lower, p_eur_upper, possible_p_eur, p_eur_pdf_vals = neanderthal_frequency_dist_with_uncertainty(p_eur, n_eur_indv)
        # p_afr_pdf, p_afr_lower, p_afr_upper, possible_p_afr, p_afr_pdf_vals = neanderthal_frequency_dist_with_uncertainty(p_afr, n_afr_indv)
        p_eas_pdf, p_eas_lower, p_eas_upper, possible_p_eas, p_eas_pdf_vals = neanderthal_frequency_dist_with_uncertainty(p_eas, n_eas_indv)
        p_amr_pdf, p_amr_lower, p_amr_upper, possible_amr, prob_amr = neanderthal_frequency_dist_with_uncertainty(n_haplotypes / n_amr_indv,
                                                                                                                  n_amr_indv, discretize=True)

        # get min and max contribution per ancestry
        # min_contr_afr, max_contr_afr = get_min_and_max_contributions(trials_afr, p_afr_lower, p_afr_upper)
        min_contr_eur, max_contr_eur = get_min_and_max_contributions(trials_eur, p_eur_lower, p_eur_upper)
        min_contr_eas, max_contr_eas = get_min_and_max_contributions(trials_eas, p_eas_lower, p_eas_upper)
        # get the unique contributions for each ancestry to avoid redundant calculations
        # unique_afr = np.arange(min_contr_afr, max_contr_afr)
        unique_eur = np.arange(min_contr_eur, max_contr_eur)
        unique_eas = np.arange(min_contr_eas, max_contr_eas)

        # calculate probabilities of unique contributions for each ancestry given all possible success rates
        prob_eur = binom.pmf(np.repeat(unique_eur, possible_p_eur.shape[0]).reshape(unique_eur.shape[0],
                                                                                    possible_p_eur.shape[0]),
                             trials_eur, possible_p_eur)
        # prob_afr = binom.pmf(np.repeat(unique_afr, possible_p_afr.shape[0]).reshape(unique_afr.shape[0],
        #                                                                             possible_p_afr.shape[0]),
        #                      trials_afr, possible_p_afr)
        prob_eas = binom.pmf(np.repeat(unique_eas, possible_p_eas.shape[0]).reshape(unique_eas.shape[0],
                                                                                    possible_p_eas.shape[0]),
                             trials_eas, possible_p_eas)

        # calculate joint probability of certain contributions and all possible success rate per ancestry
        # prob_afr = (prob_afr * p_afr_pdf_vals).sum(axis=1)
        prob_eur = (prob_eur * p_eur_pdf_vals).sum(axis=1)
        prob_eas = (prob_eas * p_eas_pdf_vals).sum(axis=1)
        # prob_afr = np.nan_to_num(prob_afr, nan=1)
        prob_eur = np.nan_to_num(prob_eur, nan=1)
        prob_eas = np.nan_to_num(prob_eas, nan=1)
        # multiply probabilities of all ancestries
        # all_prob = (prob_eur[:, np.newaxis] * prob_afr)[:, :, np.newaxis] * prob_eas
        all_prob = (prob_eas[:, np.newaxis] * prob_eur)
        # get all possible contributions
        # afr, eur, eas = np.meshgrid(unique_afr, unique_eur, unique_eas)
        eur, eas = np.meshgrid(unique_eur, unique_eas)
        expectation = (all_prob * (eur + eas)).sum()
        expectations.append(expectation)
        updated_n_haplotypes = np.rint(sum(possible_amr * prob_amr))
        corrected_n_haplotypes.append(updated_n_haplotypes)
        # don't calculate probabilities
        if not calculate_probabilities:
            continue
        # depletion
        exp_prob = pd.DataFrame.from_dict({"exp": (eur + eas).flatten(),
                                           'prob': all_prob.flatten()}).groupby('exp').sum().sort_index().cumsum()
        if expectation > possible_amr[-1]:
            comparison = 'smaller'
            thresh_amr = possible_amr[np.cumsum(prob_amr[::-1]) <= 0.975][-1]
            thresh_exp = exp_prob[exp_prob.prob > 0.025].index[0]
        # enrichment
        elif expectation < possible_amr[0]:
            comparison = 'greater'
            thresh_amr = possible_amr[np.cumsum(prob_amr[::-1]) >= 0.025][0]
            thresh_exp = exp_prob[exp_prob.prob < 0.975].index[-1]

        # ranges overlap
        elif possible_amr[0] <= expectation <= possible_amr[-1]:
            # selection_coefficients_low_ci.append(0.0)
            # selection_coefficients_high_ci.append(0.0)
            selection_coefficients.append(0.0)
            probabilities.append(1.0)
            continue
        else:
            print(c_split)
            sys.exit(1)
        # get indices of relevant combinations of contributions
        # dimensions = get_dimensions(afr, eur, eas, possible_amr[0], comparison)
        dimensions = get_dimensions(eur, eas, thresh_amr, comparison)

        # observing a value that is very unlikely given ancestry proportions and introgression frequencies in
        # reference populations
        if dimensions[0].shape[0] == 0:
            c_prob = 0
        # pick the relevant probabilities and multiply with probability of observing certain number of success
        # if lowest index is > 0 I need to renormalize to match axes of all_prob
        else:
            # c_prob = all_prob[dimensions[0],
            #                   dimensions[1],
            #                   dimensions[2]].sum()
            c_prob = all_prob[dimensions[0],
                              dimensions[1]].sum()
        selection_coefficients.append(estimate_selection_coefficient(thresh_exp / n_amr_indv, thresh_amr / n_amr_indv))
        probabilities.append(c_prob)
    if calculate_probabilities:
        df['probabilities'] = probabilities
        df['sel_coeff'] = selection_coefficients
    df['expectations'] = expectations
    df['freq_updated'] = corrected_n_haplotypes
    logging.info(f"Finished split {c_split}/{n_splits}\n")

    return df


def round_local_ancestry_haplotypes(df, n_amr_indv):
    """
    Round number of haplotypes with recent ancestries to match total number of indivduals
    :param df: pd.DataFrame
    :param n_amr_indv: int, number of individuals
    :return: pd.DataFrame
    """
    # rescale number of local ancestry haplotypes to n_amr_indv
    la_freq = df.loc[:, ['freq_la_afr', 'freq_la_eur', 'freq_la_eas']].values * \
              (n_amr_indv / df.loc[:, ['freq_la_afr', 'freq_la_eur', 'freq_la_eas']] \
               .sum(axis=1).values)[:, np.newaxis]
    # There are some edge cases when you end with n_amr_indv +/- 1 haplotype --> need to handle those
    la_freq_rounded = np.round(la_freq)
    rows_one_less = np.where(la_freq_rounded.sum(axis=1) < n_amr_indv)[0]
    rows_one_more = np.where(la_freq_rounded.sum(axis=1) > n_amr_indv)[0]

    # add one count to a random columns where we have too few haplotypes,
    # pick a column with probability equal to the rounding difference
    probs = la_freq[rows_one_less] - la_freq_rounded[rows_one_less]
    if probs.shape[0] > 0:
        cols = np.apply_along_axis(sample_cols, 1, probs)
        la_freq_rounded[rows_one_less] = np.apply_along_axis(add_count, 1,
                                                             np.hstack([la_freq_rounded[rows_one_less],
                                                                        cols]), 1)[:, :3]
    # subtract one count to a random columns where we have too many haplotypes,
    # pick a column with probability equal to the rounding difference
    probs = np.abs(la_freq[rows_one_more] - la_freq_rounded[rows_one_more])
    if probs.shape[0] > 0:
        cols = np.apply_along_axis(sample_cols, 1, probs)
        la_freq_rounded[rows_one_more] = np.apply_along_axis(add_count, 1,
                                                             np.hstack([la_freq_rounded[rows_one_more],
                                                                        cols]), -1)[:, :3]

    df.loc[:, ['freq_la_afr', 'freq_la_eur', 'freq_la_eas']] = la_freq_rounded
    return df


# def calculate_probability_of_observed_frequencies(df, n_afr_indv, n_eur_indv, n_eas_indv, n_amr_indv, threads,
#                                                   calculate_probabilities=False):
def calculate_probability_of_observed_frequencies(df, n_eur_indv, n_eas_indv, n_amr_indv, threads,
                                                  calculate_probabilities=False):
    """
    Calculate probability of observed introgressed frequency in target population given local ancestry frequencies and
    introgression frequencies in reference populations
    :param df: pd.DataFrame, Merged dataframe containing information about introgression frequency in target and
                             reference populations as well as local ancestry frequency in target population
    :param n_afr_indv: int, Number of African reference individuals
    :param n_eur_indv: int, Number of European reference individuals
    :param n_eas_indv: int, Number of East Asian reference individuals
    :param n_amr_indv: int, Number of African-American target individuals
    :param threads: int, Number of CPUs
    :param: calculate_probabilities: boolean, whether to calculate probabilities. If false just calculate expectations
    :return: pd.DataFrame, Same as input dataframe but also including expectations of introgressed frequencies and
                           probabilities of observed introgressed frequencies
    """
    window_size = df.index.get_level_values('end')[0] - df.index.get_level_values('start')[0]
    # introgression freqs in different continental super populations
    # df.freq_afr /= n_afr_indv * window_size
    df.freq_eur /= n_eur_indv * window_size
    df.freq_eas /= n_eas_indv * window_size
    df.freq = np.ceil(df.freq.values / window_size)
    #     df.freq = np.floor(df.freq.values / window_size)

    # rescale number of local ancestry haplotypes to n_amr_indv
    df = round_local_ancestry_haplotypes(df, n_amr_indv)
    df.reset_index(inplace=True)
    splits = np.array_split(df.sample(frac=1, random_state=42), df.shape[0] // 4)
    pool = mp.Pool(processes=threads)
    ready_to_map = [(split, n_eur_indv,# n_afr_indv,
                     n_eas_indv, n_amr_indv, calculate_probabilities, i, len(splits))
                    for i, split in enumerate(splits)]
    calculate_expectation_and_probability_helper(ready_to_map[0])
    dfs = pool.map(calculate_expectation_and_probability_helper, ready_to_map)
    pool.close()
    pool.join()
    df = pd.concat(dfs)
    df.sort_values(['chrom', 'start'], inplace=True)
    return df


def sample_control_sequences_helper(args):
    """
    Helper function to sample control segments in parallel
    :rtype: pd.DataFrame, control segments matched for length (and introgression frequency in EUR)
    """
    reference_df, df_to_sample_from, genome_file, seed = args
    sampled_control_segments = sample_control_sequences(reference_df, df_to_sample_from, genome_file, seed)
    return sampled_control_segments


def sample_control_sequences(reference_df, df_to_sample_from, genome_file, seed):
    """
    Sample control introgressed sequences with a similar length distribution as putatively selected introgressed
    segments from a set of introgressed segments that show no evidence of selection. Also matches for introgression
    frequency in Europe if the column freq_eur_bin is included in reference_df and df_to_sample_from.
    :param reference_df: pd.DataFrame, set of putatively selected introgressed segments
    :param df_to_sample_from: pd.DataFrame, set of introgressed segments without evidence of selections
    :param genome_file: str, path to genome file that can be used by BedTools
    :param seed: int, random seed (different seeds are needed when multiprocessing)
    :return: pd.DataFrame, control segments without evidence of selection and a similar length distribution
    """
    np.random.seed(seed)
    columns = df_to_sample_from.columns.tolist()
    control_segments = BedTool.from_dataframe(df_to_sample_from).closest(
        BedTool.from_dataframe(reference_df.sort_values(['chrom', 'start', 'end'])),
        d=True, g=genome_file).to_dataframe(names=np.concatenate([df_to_sample_from.columns.values,
                                                                  reference_df.columns.values + '_b', ['distance']]))
    # select control segments on same chromosome that are at least 0.5Mb away
    control_segments = control_segments[(control_segments.distance != -1) & (control_segments.distance >= 5e5)]
    control_seq = []
    # random_sample_length = []
    if 'freq_eur_bin' in columns:
        matches = reference_df.groupby(['length', 'freq_eur_bin'])['start'].count().reset_index().values
    else:
        matches = reference_df.groupby(['length'])['start'].count().reset_index().values
    for to_match in matches:

        # need to match for introgression frequency in Europe
        if to_match.shape[0] == 3:
            length, freq_bin, count = to_match
            min_freq_bin = freq_bin
            max_freq_bin = freq_bin
        else:
            length, count = to_match
        while count > 0:
            # if control_segments[(control_segments.chrom == chrom)].shape[0] > 0:
            if control_segments.shape[0] > 0:
                if to_match.shape[0] == 3:
                    min_freq_bin = freq_bin
                    max_freq_bin = freq_bin
                min_length = length
                # max_length = length
                control = None
                # try sampling a segment with exact same length, otherwise give a length range
                while control is None:
                    try:
                        # need to match for introgression frequency in Europe
                        if to_match.shape[0] == 3:
                            control = control_segments[(control_segments.length >= min_length) &
                                                       (control_segments.freq_eur_bin >= min_freq_bin) &
                                                       (control_segments.freq_eur_bin <= max_freq_bin)].sample(n=1,
                                                                                                               replace=False)

                        else:
                            control = control_segments[(control_segments.length >= min_length)].sample(n=1,
                                                                                                       replace=False)
                        start = np.random.randint(control.start, control.end - length - 1)
                        end = np.random.randint(start + min_length, start + length + 1)
                        control.start = start
                        control.end = end
                    except ValueError:
                        # ensure there are always some control segments fulfilling length requirement
                        if to_match.shape[0] == 2 and control_segments[(control_segments.length >= min_length)].shape[0] == 0:
                            min_length -= 5000
                        # max_length += 5000
                        # ensure there are always some control segments fulfilling freq bin requirement
                        if to_match.shape[0] == 3 and control_segments[(control_segments.length >= min_length) &
                                                                       (control_segments.freq_eur_bin >= min_freq_bin) &
                                                                       (control_segments.freq_eur_bin <= max_freq_bin)].shape[0] == 0:
                            # prioritize freq bin over length
                            if control_segments[(control_segments.freq_eur_bin >= min_freq_bin) &
                                                (control_segments.freq_eur_bin <= max_freq_bin)].shape[0] == 0:
                                min_freq_bin -= 1
                                max_freq_bin += 1
                                # reset length
                                min_length = length
                            else:
                                min_length -= 5000
                    control_seq.append(control)
                # drop sampled segment and segments possibly in LD to avoid biases
                control_segments = control_segments[~((control_segments.chrom == control.chrom.values[0]) &
                                                      (control_segments.end > control.start.values[0] - 500000) &
                                                      (control_segments.start < control.end.values[0] + 500000))]
                count -= 1
            else:
                break
    try:
        sampled_control_segments = pd.concat(control_seq)
    except ValueError:
        sampled_control_segments = pd.DataFrame(columns=columns)
    if 'freq_eur_bin' in columns:
        columns.remove('freq_eur_bin')
    sampled_control_segments = sampled_control_segments.loc[:, columns].sort_values(['chrom', 'start', 'end'])

    return sampled_control_segments


def estimate_effective_number_of_independent_tests_helper(args):
    """
    Helper function to estimate effective number of independent tests chromosome-wise
    :rtype: float, effective number of tests
    """
    df_coverage, c_df_freq, aa_individuals, n_amr_indv, chrom_size, max_eas_anc, min_expectation, window_size, step_size, method = args
    # get windows_to_include
    windows_to_include = np.where((c_df_freq.freq_afr == 0).values &
                                  (c_df_freq.freq_la_eas < np.floor(max_eas_anc * n_amr_indv)).values &
                                  (c_df_freq.freq_eur > 0).values &
                                  (c_df_freq.expectations >= min_expectation).values)[0]
    # calculate window midpoints
    window_midpoints = np.arange(0, chrom_size + step_size, step_size) + window_size / 2
    # initialize coverage matrix
    coverage_mat = np.zeros((len(aa_individuals), len(window_midpoints)))
    # iterate over all african american individuals
    for i, indv in enumerate(aa_individuals):
        # pick them and corresponding chromosome
        c_df = df_coverage[df_coverage.iid == indv]
        if c_df.shape[0] == 0:
            continue
        # iterate over all windows at which they show neanderthal introgression and save coverage
        windows = (c_df.start.values + c_df.end.values) / 2
        for win, cov in zip(windows, c_df.coverage.values):
            try:
                j = np.where(window_midpoints == win)[0][0]
            except IndexError:
                # some edge cases at the end of chromosome
                j = np.where(np.isclose(window_midpoints, win, rtol=0, atol=window_size / 2))[0][0]
            coverage_mat[i, j] = cov
    # calculate physical distance between windows
    distance_matrix = pairwise_distances(window_midpoints.reshape(-1, 1)) / 1000
    distance_matrix = distance_matrix[:, coverage_mat.sum(axis=0) > 0][coverage_mat.sum(axis=0) > 0, :]
    distance_matrix = distance_matrix[:, windows_to_include][windows_to_include, :]
    # only pick windows that show evidence of introgression in any window
    coverage_mat = coverage_mat[:, coverage_mat.sum(axis=0) > 0]
    # only pick windows that pass other filtering (see docstring)
    coverage_mat = coverage_mat[:, windows_to_include]
    if method == 'all':
        # assume all windows are independent
        meff = coverage_mat.shape[1]
    elif method == 'gao' or method == 'li_and_ji' or method == 'nyholt' or method == 'galeway':
        # calculate correlation matrix of coverage across windows
        correlation_matrix = pd.DataFrame(coverage_mat).corr()
        # set correlation to 0 for windows separated by more than 1 Mb
        correlation_matrix = np.where(distance_matrix > 1000, 0, correlation_matrix)
        # calculate eigenvalues
        eigenvalues = np.real(np.linalg.eigvals(correlation_matrix))
        # estimate effective number of independent tests using specified method
        if method == 'gao':
            meff = meff_gao(eigenvalues)
        elif method == 'li_and_ji':
            meff = meff_li_ji(eigenvalues)
        elif method == 'nyholt':
            meff = meff_nyholt(eigenvalues)
        elif method == 'galeway':
            meff = meff_galeway(eigenvalues)
    else:
        raise NotImplementedError(f"{method} for estimating the effective of independent tests is not supported")
    return meff


def estimate_effective_number_of_independent_tests(df_coverage, df_freq, aa_individuals, n_amr_indvs, chrom_sizes,
                                                   max_eas_anc, min_expectation, window_size, step_size, method, threads):
    """
    Estimate effective number of independent tests using specified method
    :param df_coverage: pd.DataFrame, introgression coverage per window per individual
    :param df_freq: pd.DataFrame, introgression frequency per region and other statistics
    :param aa_individuals: list-like, IIDs of African-American windows
    :param n_amr_indvs: int, number of African American genomes
    :param chrom_sizes: dict, chromosome sizes
    :param max_eas_anc: float, maximal EAS ancestry
    :param min_expectation: int, minimum number of expected introgressed haplotypes overlapping a given window
    :param window_size: int, window size used to tile genome
    :param step_size: int, step size used to tile genome
    :param method: str, method to be used for estimating effective number of independent tests.
                        One of all, gao, li_and_ji, nyholt, or galeway
    :param threads: int, number of CPUs
    :return: int, floor of number of estimated effective independent tests
    """
    ready_to_map = [(df_coverage[df_coverage.chrom == chrom], df_freq[df_freq.chrom == chrom], aa_individuals,
                     n_amr_indvs, chrom_sizes[chrom], max_eas_anc, min_expectation, window_size, step_size, method)
                    for chrom in df_coverage.chrom.unique()]
    pool = mp.Pool(processes=threads)
    meff = pool.map(estimate_effective_number_of_independent_tests_helper, ready_to_map)
    pool.close()
    pool.join()
    return np.floor(sum(meff))


def calculate_combined_pvalues(df, original_df, method='harmonic'):
    """
    Combine p-values according to Fisher's method or use harmonic mean
    :param df: pd.DataFrame, merged segments for which to calculate combined p-values from original_df
    :param original_df: pd.DataFrame, original dataframe with one (FDR-corrrected) p-value for each window
    :param method: str, specifing according to which method to combine p-values
    :return: pd.DataFrame, df with extra columns of combined (FDR-corrected) p-values for each segment
    """
    probabilities = []
    corrected = []
    for chrom, start, end in zip(df.chrom, df.start, df.end):
        c_df = original_df[(original_df.chrom == chrom) & (original_df.start >= start) & (original_df.end <= end)]
        if method == 'fisher':
            probabilities.append(combine_pvalues(c_df.probabilities.values).pvalue)
            corrected.append(combine_pvalues(c_df.corrected.values).pvalue)
        elif method == 'harmonic':
            probabilities.append(hmean(c_df.probabilities.values))
            corrected.append(hmean(c_df.corrected.values))
        else:
            raise NotImplementedError
    df['probabilities'] = probabilities
    df['corrected'] = corrected
    return df


def subtract_mask(args):
    """
    Helper function to subtract mask in parallele. args should be of form (segments, masks)
    :return: pd.DataFrame, masked segments
    """
    df = args[0].copy()
    df.chrom = [int(chrom.replace('chr', '')) for chrom in df.chrom.values]
    df_masked = BedTool.from_dataframe(df).subtract(args[1]).to_dataframe()
    df_masked.chrom = [f'chr{chrom}' for chrom in df_masked.chrom.values]
    return df_masked


def join_data_frames_by_chromosomes(args):
    """
    Function to join different dataframes chromosome-wise in parallel
    :param args: tuple, nea_amr_freq, nea_afr_freq, nea_eur_freq, nea_eas_freq, local_anc_freq
    :rtype: pd.DataFrame, joined dataframe
    """
    # nea_amr_freq, nea_afr_freq, nea_eur_freq, nea_eas_freq, local_anc_freq = args
    nea_amr_freq, nea_eur_freq, nea_eas_freq, local_anc_freq = args

    # df = nea_amr_freq.set_index(['chrom', 'start', 'end']).join(nea_afr_freq.set_index(['chrom', 'start', 'end']),
    #                                                             rsuffix='_afr')
    df = nea_amr_freq.set_index(['chrom', 'start', 'end']).join(nea_eur_freq.set_index(['chrom', 'start', 'end']),
                                                                rsuffix='_eur', how='outer')

    # df = df.join(nea_eur_freq.set_index(['chrom', 'start', 'end']), rsuffix='_eur')
    df = df.join(nea_eas_freq.set_index(['chrom', 'start', 'end']), rsuffix='_eas', how='outer')
    df.fillna(0.0, inplace=True)
    df = df.join(local_anc_freq.loc[local_anc_freq.la == 'AFR', ['chrom', 'start', 'end', 'freq']]. \
                 set_index(['chrom', 'start', 'end']), rsuffix='_la_afr')
    df = df.join(local_anc_freq.loc[local_anc_freq.la == 'EUR', ['chrom', 'start', 'end', 'freq']]. \
                 set_index(['chrom', 'start', 'end']), rsuffix='_la_eur')
    df = df.join(local_anc_freq.loc[local_anc_freq.la == 'EAS', ['chrom', 'start', 'end', 'freq']]. \
                 set_index(['chrom', 'start', 'end']), rsuffix='_la_eas')
    return df


def get_recombination_mask(genetic_maps_fn, chrom_sizes_fn, genome_file_fn):
    """
    Get windows with "average" recombination rate --> middle 33%
    :param genetic_maps_fn: str, paths to genetic maps for all chromosomes
    :param chrom_sizes_fn: str, path to file with chromosome sizes
    :param genome_file_fn: str, path to genome file as bedtools expects it
    :return: df, windows that pass recombination filter
    """
    genetic_maps = []
    # load genetic maps
    for gm in genetic_maps_fn:
        genetic_map = pd.read_csv(gm, sep='\t')
        genetic_map.rename({'Chromosome': 'chrom_map', 'Position(bp)': 'end_map', 'Rate(cM/Mb)': 'rate'}, axis=1,
                           inplace=True)
        genetic_map['start_map'] = np.concatenate([[0], genetic_map.end_map[:-1]])
        genetic_map = genetic_map.loc[:, ['chrom_map', 'start_map', 'end_map', 'rate']]
        genetic_maps.append(genetic_map)
    genetic_map = pd.concat(genetic_maps)
    # window genome into windows of 300kb
    windows = BedTool(chrom_sizes_fn).window_maker(w=300000, g=genome_file_fn)
    columns = ['chrom', 'start', 'end']
    columns.extend(genetic_map.columns.values.tolist())
    columns.append('overlap')
    # calculate rec rate per window
    rec_windows = windows.intersect(BedTool.from_dataframe(genetic_map), wao=True).to_dataframe(names=columns)
    rec_windows.rate = np.where(rec_windows.rate == '.', '0', rec_windows.rate)
    rec_windows.rate = rec_windows.rate.astype(float)
    rec_windows['rate'] *= rec_windows['overlap']
    rec_windows = rec_windows.loc[:, ['chrom', 'start', 'end',
                                      'rate', 'overlap']].groupby(['chrom', 'start', 'end']).sum()
    rec_windows.rate /= rec_windows.overlap
    rec_windows.reset_index(inplace=True)
    rec_windows.fillna(0, inplace=True)
    # select windows in the middle 33%
    low_rec_thresh = np.percentile(rec_windows.rate.values, 33)
    high_rec_thresh = np.percentile(rec_windows.rate.values, 66)
    pass_rec_filter = rec_windows.loc[(rec_windows.rate >= low_rec_thresh) &
                                      (rec_windows.rate <= high_rec_thresh), ['chrom', 'start', 'end']]
    return pass_rec_filter


def find_significant_segments(df, enriched=False, depleted=False, min_length=90000, windowsize=50000, stride=10000):
    if enriched:
        df = df[df.freq > df.expectations].copy()
    elif depleted:
        df = df[df.freq < df.expectations].copy()
    df.reset_index(drop=True, inplace=True)
    if windowsize > stride:
        min_windows = int((min_length - windowsize) / stride + 1)
    elif windowsize == stride:
        min_windows = int(min_length / stride)
    idx_start = 0
    idx_end = 1
    differences = np.diff(df.start.values)
    to_keep = []
    for diff in differences:
        if diff == stride:
            idx_end += 1
        elif diff != stride and idx_end - idx_start >= min_windows:
            to_keep.append(np.arange(idx_start, idx_end + 1))
            idx_start = idx_end
            idx_end += 1
        else:
            idx_start = idx_end
            idx_end += 1
    if idx_end - idx_start >= min_windows:
        to_keep.append(np.arange(idx_start, idx_end))
    if len(to_keep) > 0 and df.shape[0] > 0:
        to_keep = np.concatenate(to_keep)
        # need to trait pvalues differently
        columns = np.where((df.columns.values != 'probabilities') & (df.columns.values != 'corrected'))[0][3:] + 1
        merged = BedTool.from_dataframe(df.iloc[to_keep]).merge(c=columns.tolist(),
                                                                o='mean',
                                                                header=True)
        df.drop(['probabilities', 'corrected'], axis=1, inplace=True)
        merged = merged.to_dataframe(names=df.columns)
    else:
        merged = pd.DataFrame(columns=df.columns)
    return merged


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--lai', nargs='+', help='Local ancestry information inferred using flare for '
                                                       'all chromosomes')
    parser.add_argument('--nea_amr', help='Introgression frequencies in African-Americans')
    #parser.add_argument('--nea_afr', help='Introgression frequencies in Africans')
    parser.add_argument('--nea_eur', help='Introgression frequencies in Europeans')
    parser.add_argument('--nea_eas', help='Introgression frequencies in East Asians')
    parser.add_argument('--ibdmix_all', help='IBDmix results for all populations')
    parser.add_argument('--ibdmix_all_coverage', help='IBDmix results with introgression frequencies for '
                                                      'each window and individual')
    parser.add_argument('--max_eas_anc', type=float, help='Maximal East Asian ancestry at a given position')
    parser.add_argument('--min_length', type=int, help='Minimal length of introgressed segments')
    parser.add_argument('--min_expectation', type=int, help='Minimal expected number of introgressed haplotypes '
                                                            'overlapping a window')
    parser.add_argument('-g', '--genome_file', help="'Genome file' for bedtools indicating chromosome sizes")
    parser.add_argument('--windowed_genome_file', help="'Windowed genome file' for bedtools indicating chromosome sizes")
    parser.add_argument('--genetic_maps', nargs='+', help='Paths to genetic maps in plink format')
    parser.add_argument('--chrom_sizes', help='Path to file with chromosome sizes')
    parser.add_argument('-os', '--output_selected', help='Output filename which will contain putatively selected '
                                                         'Neanderthal segments in African-Americans')
    parser.add_argument('-ons', '--output_not_selected', help='Prefix for output filen that will contain not selected '
                                                              'control Neanderthal segments in African-Americans. '
                                                              '"_{n}.bed" suffix will be appended')
    parser.add_argument('-oip', '--output_ibdmix_all_pvalues', help='Output file name for file that will be same as '
                                                                    '--ibdmix_all_coverage but will also contain '
                                                                    'corrected pvalues for windows of interest')
    parser.add_argument('-oie', '--output_ibdmix_all_expectations', help='Output file name for file that will be same '
                                                                         'as --ibdmix_all_coverage but will also '
                                                                         'contain expected introgression frequencies')
    parser.add_argument('-w', '--window_size', type=int, help='Window size of bp used to calculate test statistic')
    parser.add_argument('-s', '--step_size', type=int, help='Step size of bp used to calculate test statistic')

    parser.add_argument('--meff_method', help='Method to estimate effective number of independent tests. One of all, '
                                              'gao, li_and_ji, nyholt, or galeway. default=li_and_ji',
                        default='li_and_ji')
    parser.add_argument('--alpha', type=float, help='FDR error rate')
    parser.add_argument('--reps', type=int, default=100,
                        help='Number of bootstrap sets to create from control segments. [100]')
    parser.add_argument('--n_amr', type=int, help='Number of admixed African American individuals. '
                                                  'If not set will be inferred from data. [None]',
                        default=None)
    parser.add_argument('--n_eas', type=int, help='Number of East Asian individuals. '
                                                  'If not set will be inferred from data. [None]',
                        default=None)
    parser.add_argument('--n_afr', type=int, help='Number of African individuals. '
                                                  'If not set will be inferred from data. [None]',
                        default=None)
    parser.add_argument('--n_eur', type=int, help='Number of European individuals. '
                                                  'If not set will be inferred from data. [None]',
                        default=None),
    parser.add_argument('--chrom', help='Consider only a specific chromosome. [None]', default=None)
    parser.add_argument('-t', '--threads', help='Number of CPUs. [1]', type=int)
    parser.add_argument('--multiple_testing_correction', help='Which method to use for multiple testing '
                                                              '(fdr or bonferroni). bonferroni',
                        default='bonferroni')
    parser.add_argument('--mask_pattern', help='String pattern of mask files, must allow string formatting of "chrom"')
    parser.add_argument('--tmp_dir', help='Temporary directory for pybedtools. Default=/tmp/', default='/tmp/')

    args = parser.parse_args()
    logging.basicConfig(filename='identify_putatively_selected_introgressed_segments.log', level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s', datefmt='%H:%M:%S')

    logger = logging.getLogger()
    install_mp_handler()
    lai = args.lai
    nea_amr = args.nea_amr
    #nea_afr = args.nea_afr
    nea_eur = args.nea_eur
    nea_eas = args.nea_eas
    ibdmix_all = args.ibdmix_all
    ibdmix_all_coverage = args.ibdmix_all_coverage

    # read data
    logging.info('Loading local ancestry data')
    local_anc_freq = pd.concat([pd.read_csv(flare, sep='\t', names=['chrom', 'start', 'end', 'la', 'freq'])
                                for flare in lai])
    logging.info('Loading introgression frequencies data')
    #nea_afr_freq = pd.read_csv(nea_afr, sep='\t', names=['chrom', 'start', 'end', 'freq'], usecols=[0, 1, 2, 4])
    nea_amr_freq = pd.read_csv(nea_amr, sep='\t', names=['chrom', 'start', 'end', 'freq'], usecols=[0, 1, 2, 4])
    nea_eur_freq = pd.read_csv(nea_eur, sep='\t', names=['chrom', 'start', 'end', 'freq'], usecols=[0, 1, 2, 4])
    nea_eas_freq = pd.read_csv(nea_eas, sep='\t', names=['chrom', 'start', 'end', 'freq'], usecols=[0, 1, 2, 4])
    ibdmix = pd.read_csv(ibdmix_all, sep='\t', names=['chrom', 'start', 'end', "LOD", "IID", "pop", "super_pop"])
    df_coverage = pd.read_csv(ibdmix_all_coverage, sep='\t', names=['chrom', 'start', 'end', 'iid', 'coverage'])
    df_coverage = df_coverage.join(ibdmix.loc[:, ['IID', 'super_pop']].drop_duplicates().set_index('IID'), on='iid')
    # in simulations we only have AA
    if not 'AMR' in ibdmix.super_pop.unique():
        df_coverage = df_coverage[df_coverage.super_pop == 'AA']
    else:
        df_coverage = df_coverage[df_coverage.super_pop == 'AMR']
    df_coverage.drop('super_pop', axis=1, inplace=True)
    # join dataframes
    logging.info('Joining local ancestry data and introgression frequencies')
    read_to_map_join_df = [(nea_amr_freq[nea_amr_freq.chrom == chrom],# nea_afr_freq[nea_afr_freq.chrom == chrom],
                            nea_eur_freq[nea_eur_freq.chrom == chrom], nea_eas_freq[nea_eas_freq.chrom == chrom],
                            local_anc_freq[local_anc_freq.chrom == chrom]) for chrom in local_anc_freq.chrom.unique()]
    pool = mp.Pool(processes=args.threads)
    dfs = pool.map(join_data_frames_by_chromosomes, read_to_map_join_df)
    pool.close()
    pool.join()
    df = pd.concat(dfs)
    df.sort_values(['chrom', 'start', 'end'], inplace=True)
    df.fillna(0.0, inplace=True)

    # if args.n_afr is None:
    #     n_afr_indv = ibdmix.loc[ibdmix.super_pop == 'AFR', 'IID'].unique().shape[0]
    # else:
    #     n_afr_indv = args.n_afr
    if args.n_amr is None:
        if not 'AMR' in ibdmix.super_pop.unique():
            n_amr_indv = ibdmix.loc[ibdmix.super_pop == 'AA', 'IID'].unique().shape[0]
        else:
            n_amr_indv = ibdmix.loc[ibdmix.super_pop == 'AMR', 'IID'].unique().shape[0]
    else:
        n_amr_indv = args.n_amr
    if args.n_eur is None:
        n_eur_indv = ibdmix.loc[ibdmix.super_pop == 'EUR', 'IID'].unique().shape[0]
    else:
        n_eur_indv = args.n_eur
    if args.n_eas is None:
        n_eas_indv = ibdmix.loc[ibdmix.super_pop == 'EAS', 'IID'].unique().shape[0]
    else:
        n_eas_indv = args.n_eas
    df_subsetted = df.copy()
    logging.info('Calculating expected introgression_frequencies for all windows')
    # Calculate expectations only
    df = calculate_probability_of_observed_frequencies(df,# n_afr_indv,
                                                       n_eur_indv, n_eas_indv, n_amr_indv, args.threads)
    # save to dataframe
    # df.chrom = df.chrom.astype(int)
    df.start = df.start.astype(int)
    df.end = df.end.astype(int)
    df.to_csv(args.output_ibdmix_all_expectations, index=False, header=True, sep='\t')
    # Subset to regions of interest for computational efficiency --> don't need to do pvalue calculation for all windows
    # Also reduces multiple testing burden
    logging.info('Calculating pvalues for introgression_frequencies for subsetted windows windows')

    df_subsetted.reset_index(inplace=True)

    windows = pd.read_csv(args.windowed_genome_file, sep='\t', names=['chrom', 'start', 'end'])
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
    ready_to_map_mask = [(windows, mask, chrom) for chrom in windows.chrom.unique()]
    pool = mp.Pool(processes=args.threads)
    masked = pool.map(get_overlap_mask, ready_to_map_mask)
    pool.close()
    pool.join()
    masked = pd.concat(masked)
    masked.sort_values(['chrom', 'start'])
    columns = df.columns.values.tolist()
    columns.extend([col + "_m" for col in masked.columns.values])
    columns.append("overlap")
    df = BedTool.from_dataframe(df).intersect(BedTool.from_dataframe(masked),
                                              wo=True, f=1, r=True).to_dataframe(names=columns)
    df.drop(['chrom_m', 'start_m', "end_m", "overlap"], axis=1, inplace=True)
    df.rename({'masked_m': "masked"}, axis=1, inplace=True)
    df.masked /= args.window_size

    logging.info('Calculating pvalues for introgression_frequencies for subsetted windows windows')

    df_subsetted = df_subsetted[(df.expectations >= args.min_expectation) &
                                # (df.freq_afr == 0) &
                                (df.freq_la_eur >= np.ceil(0.1 * n_amr_indv)) &
                                (df.freq_la_afr >= np.ceil(0.5 * n_amr_indv)) &
                                (df.freq_la_eas < np.floor(args.max_eas_anc * n_amr_indv)) &
                                # (df.freq_eur > 0) &
                                (df.masked < 0.5)]

    df_subsetted.set_index(['chrom', 'start', 'end'], inplace=True)
    df_subsetted = calculate_probability_of_observed_frequencies(df_subsetted,# n_afr_indv,
                                                                 n_eur_indv, n_eas_indv,
                                                                 n_amr_indv, args.threads, calculate_probabilities=True)
    # # Apply recombination mask
    pass_rec_filter = get_recombination_mask(args.genetic_maps, args.chrom_sizes, args.genome_file)
    df_subsetted = BedTool.from_dataframe(df_subsetted).intersect(
        BedTool.from_dataframe(pass_rec_filter), f=1).to_dataframe(names=df_subsetted.columns)
    # Calculate expected variance from genetic drift over the last 15 generations
    sd_drift = lambda x, n: np.sqrt((x / n) * (1 - (x / n)) * (1 - np.exp(-15 / 30000)))
    lower_bound = t.ppf(0.025 / df_subsetted.shape[0], n_amr_indv,
                        loc=df_subsetted.expectations.values / n_amr_indv,
                        scale=sd_drift(df_subsetted.expectations.values, n_amr_indv))
    lower_bound = np.nan_to_num(np.where(lower_bound < 0, 0, lower_bound))
    upper_bound = t.ppf(1 - 0.025 / df_subsetted.shape[0], n_amr_indv,
                        loc=df_subsetted.expectations.values / n_amr_indv,
                        scale=sd_drift(df_subsetted.expectations.values, n_amr_indv))
    upper_bound = np.nan_to_num(upper_bound)
    df_subsetted['drift_lower'] = lower_bound
    df_subsetted['drift_upper'] = upper_bound
    # Correct for multiple testing
    if args.multiple_testing_correction == 'fdr':
        logging.info('Performing FDR multiple testing correction')

        # calculate effective number of independent tests for multiple testing correction
        chrom_sizes = parse_genome_file(args.genome_file)
        if not str(df.chrom.unique()[0]).startswith('chr'):
            chrom_sizes = {int(key.replace('chr', '')): val for key, val in chrom_sizes.items()}

        if not 'AMR' in ibdmix.super_pop.unique():
            meff = estimate_effective_number_of_independent_tests(df_coverage, df_subsetted,
                                                                  ibdmix.loc[ibdmix.super_pop == 'AA', 'IID'].unique(),
                                                                  n_amr_indv, chrom_sizes, args.max_eas_anc,
                                                                  args.min_expectation, args.window_size, args.step_size,
                                                                  args.meff_method, args.threads)
        else:
            meff = estimate_effective_number_of_independent_tests(df_coverage, df_subsetted,
                                                                  ibdmix.loc[ibdmix.super_pop == 'AMR', 'IID'].unique(),
                                                                  n_amr_indv, chrom_sizes, args.max_eas_anc,
                                                                  args.min_expectation, args.window_size, args.step_size,
                                                                  args.meff_method, args.threads)
        # enforce minimum number of tests, i.e., at least one test for every 150kb
        genome_size = sum([chrom_sizes[chrom] for chrom in df.chrom.unique()])
        if meff < genome_size / 150000:
            meff = np.floor(genome_size / 150000)
        print(f"MEFF:{ meff}")
        # perform FDR correction
        df_subsetted['corrected'] = fdrcorrection_meff(df_subsetted.probabilities, meff, alpha=args.alpha)[1]
    elif args.multiple_testing_correction == 'bonferroni':
        logging.info('Performing Bonferroni multiple testing correction')

        # bonferroni correction
        df_subsetted['corrected'] = np.where(df_subsetted.probabilities.values * df_subsetted.shape[0] > 1.0, 1.0,
                                             df_subsetted.probabilities.values * df_subsetted.shape[0])

    logging.info('Selecting selected regions and calculate merged pvalues')
    selected = df_subsetted[(df_subsetted.corrected < args.alpha) &
                            ((df_subsetted.freq_updated / n_amr_indv >
                              df_subsetted.drift_upper) |
                             (df_subsetted.freq_updated / n_amr_indv <
                              df_subsetted.drift_lower))]
    # select not selected segments
    not_selected = df_subsetted[(df_subsetted.corrected >= args.alpha * 2) |
                                ((df_subsetted.freq_updated / n_amr_indv >=
                                  df_subsetted.drift_lower) &
                                 (df_subsetted.freq_updated / n_amr_indv <=
                                  df_subsetted.drift_upper))]
    # merge overlapping intervals
    if selected.shape[0] == 0:
        print('No evidence of secondary selection of Neanderthal DNA')
        selected['length'] = []
        selected.to_csv(args.output_selected, index=False, header=True, sep='\t')
        if args.reps > 1:
            for i in range(args.reps):
                selected.to_csv(f'{args.output_not_selected}_{i}.bed', index=False, header=True, sep='\t')
        elif args.reps == 1:
            selected.to_csv(f'{args.output_not_selected}.bed', index=False, header=True, sep='\t')
        df_subsetted.start = df_subsetted.start.astype(int)
        df_subsetted.end = df_subsetted.end.astype(int)
        df_subsetted.to_csv(args.output_ibdmix_all_pvalues, index=False, header=True, sep='\t')
        sys.exit(0)
    df_subsetted.start = df_subsetted.start.astype(int)
    df_subsetted.end = df_subsetted.end.astype(int)
    df_subsetted.to_csv(args.output_ibdmix_all_pvalues, index=False, header=True, sep='\t')

    selected_enriched_merged = find_significant_segments(selected.sort_values(['chrom', 'start']), enriched=True,
                                                         min_length=args.min_length,  windowsize=args.window_size,
                                                         stride=args.step_size)
    selected_depleted_merged = find_significant_segments(selected.sort_values(['chrom', 'start']), depleted=True,
                                                         min_length=args.min_length, windowsize=args.window_size,
                                                         stride=args.step_size)
    selected_merged = pd.concat([selected_enriched_merged, selected_depleted_merged]).sort_values(['chrom', 'start'])

    # Fisher's method assumes independence though --> harmonic mean might be more appropriate
    # although it is still anti-conservative
    selected_merged = calculate_combined_pvalues(selected_merged, df_subsetted)
    selected_merged['length'] = selected_merged.end - selected_merged.start

    not_selected_merged = find_significant_segments(not_selected.sort_values(['chrom', 'start']),
                                                    min_length=args.min_length, windowsize=args.window_size,
                                                    stride=args.step_size)
    not_selected_merged = calculate_combined_pvalues(not_selected_merged, df_subsetted)
    not_selected_merged['length'] = not_selected_merged.end - not_selected_merged.start

    # selected_merged.chrom = selected_merged.chrom.astype(int)
    selected_merged.start = selected_merged.start.astype(int)
    selected_merged.end = selected_merged.end.astype(int)

    selected_merged.to_csv(args.output_selected, index=False, header=True, sep='\t')
    if selected_merged.shape[0] == 0:
        print('No evidence of secondary selection of Neanderthal DNA')
        if args.reps > 1:
            for i in range(args.reps):
                selected_merged.to_csv(f'{args.output_not_selected}_{i}.bed', index=False, header=True, sep='\t')
        elif args.reps == 1:
            selected_merged.to_csv(f'{args.output_not_selected}.bed', index=False, header=True, sep='\t')
        sys.exit(0)
    logging.info('Sampling non-selected control segments')

    # get bins for introgression frequency in Europe
    freq_bins = [np.percentile(selected_merged.freq_eur, x) for x in np.arange(0, 102, 2)]
    freq_bins_selected = np.zeros(selected_merged.shape[0])
    freq_bins_not_selected = np.zeros(not_selected_merged.shape[0])

    for i in range(1, len(freq_bins)):
        freq_bins_selected = np.where((selected_merged.freq_eur >= freq_bins[i - 1]) &
                                      (selected_merged.freq_eur < freq_bins[i]), i, freq_bins_selected)
        freq_bins_not_selected = np.where((not_selected_merged.freq_eur >= freq_bins[i - 1]) &
                                          (not_selected_merged.freq_eur < freq_bins[i]), i, freq_bins_not_selected)
    selected_merged['freq_eur_bin'] = freq_bins_selected
    not_selected_merged['freq_eur_bin'] = freq_bins_not_selected
    if args.reps > 1:
        ready_to_map_not_selected = [(selected_merged, not_selected_merged, args.genome_file,
                                      np.random.randint(0, 1e7)) for n in range(args.reps)]
        pool = mp.Pool(processes=args.threads)
        introgressed_control_segments_replicates = pool.map(sample_control_sequences_helper, ready_to_map_not_selected)
        pool.close()
        pool.join()
        # bootstrap sample control segments
        for i, introgressed_control_segments in \
                enumerate(introgressed_control_segments_replicates):
            introgressed_control_segments.start = introgressed_control_segments.start.astype(int)
            introgressed_control_segments.end = introgressed_control_segments.end.astype(int)
            introgressed_control_segments.to_csv(f'{args.output_not_selected}_{i}.bed', index=False, header=True,
                                                 sep='\t')
    elif args.reps == 1:
        introgressed_control_segments = sample_control_sequences(selected_merged, not_selected_merged, args.genome_file,
                                                                 42)
        introgressed_control_segments.start = introgressed_control_segments.start.astype(int)
        introgressed_control_segments.end = introgressed_control_segments.end.astype(int)
        introgressed_control_segments.to_csv(f'{args.output_not_selected}.bed', index=False, header=True, sep='\t')


if __name__ == '__main__':
    main(sys.argv)
