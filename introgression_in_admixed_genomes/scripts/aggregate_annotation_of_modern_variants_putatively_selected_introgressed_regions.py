#!/usr/bin/env python
import numpy as np
import pandas as pd
from collections import defaultdict
import sys
import argparse


def calculate_allele_frequencies(df, individuals):
    """
    Calculate allele frequencies in group of individuals
    :param df: pd.DataFrame, containing genotype information with individual IDs as columns
    :param individuals: list, individuals IDs for which to calculate allele frequencies
    :return: np.array, allele frequencies
    """
    gts = df.loc[:, individuals].replace({'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}).values
    af = gts.sum(axis=1) / (2 * gts.shape[1])
    return af


def concat_list_of_strings_in_aggregation(series):
    """
    Concatenate list of strings into a single string during pandas groupby operation
    :param series: pd.Series, values to aggregate
    :return: str, aggregated string
    """
    raw_strings = [val[2:-2].split("', '") for val in series]
    unique_strings = set(item for sublist in raw_strings for item in sublist)
    return ','.join(unique_strings)


def concat_strings_in_aggregation(series):
    """
    Concatenate strings into a single string during pandas groupby operation
    :param series: pd.Series, values to aggregate
    :return: str, aggregated string
    """
    raw_strings = [val for val in series if val != '-']
    unique_strings = set(raw_strings)
    return ','.join(unique_strings)


def concat_strings_in_aggregation_all(series):
    """
    Concatenate strings into a single string during pandas groupby operation but keep all values
    --> that is when non-repetitive
    :param series: pd.Series, values to aggregate
    :return: str, aggregated string
    """
    unique_strings = set(val for val in series if val != '-')
    return ','.join(unique_strings)


def onehotencoding_variant_annotations(series):
    """
    Perform one hot encoding of variant class annotations, e.g., VEP consequences, CLINVAR molecular consequences etc.
    :param series: pd.Series, values to encode
    :return: list, np.array, list of used categories, one hot encoded values
    """
    # Gather unique categories excluding certain values
    unique_categories = set(category for annotation in series for category in annotation
                            if category not in {'-', 'not_provided', 'not_specified'})

    # Ensure 'unknown' category is included
    unique_categories.add('unknown')

    categories = list(unique_categories)

    # Create a dictionary to map categories to indices
    onehotenc = {cat: i for i, cat in enumerate(categories)}
    encoding = []
    for annotation in series:
        cat_counts = np.zeros(len(categories))
        for csq in annotation:
            if csq in onehotenc:
                cat_counts[onehotenc[csq]] += 1
            else:
                cat_counts[onehotenc['unknown']] += 1

        encoding.append(cat_counts)
    return categories, np.vstack(encoding)


def get_basic_field_annotation(info, field, missing_value, return_int=False, return_list=False):
    """
    Helper function to easily extract infos from specific annotation field. Allow to specify return type
    and value if field is not present
    :param info: str, INFO field from VCF
    :param field: str, target annotation field in info
    :param missing_value: str, list, int, missing value to be returned when field is not present
    :param return_int: boolean, return integer instead of string
    :param return_list: boolean, return list instead of string
    :return: str, list, int, parsed annotation field value or specified missing value
    """
    try:
        if not return_int and not return_list:
            return info.split(f'{field}=')[1].split(';')[0]
        elif return_int:
            return int(info.split(f'{field}=')[1].split(';')[0])
        elif return_list:
            return [info.split(f'{field}=')[1].split(';')[0]]
    except IndexError:
        return missing_value


def parse_modern_variant_annotations(df, archaic_genomes):
    """
    Parse modern variant annotations using the fields specified below
    :param df: pd.DataFrame, BED file of putatively selected introgressed segments annotated with basic stats
    and extended by modern variant annotations
    :param archaic_genomes: list, archaic genomes that were considered during analysis
    :return: pd.DataFrame, same dataframe as input but modern variant annotations are aggregated per
                           putatively selected introgressed segment
    """
    fields = defaultdict(list)
    for ref, alt, info in zip(df.REF.values, df.ALT.values, df.INFO.values):
        try:
            ancestral = info.split('ancestral=')[1].split(';')[0]
            fields['ancestral'].append(info.split('ancestral=')[1].split(';')[0])
        except IndexError:
            ancestral = None
            fields['ancestral'].append('-')
        for archaic_genome in archaic_genomes:
            try:
                ref = info.split(f'{archaic_genome}_REF=')[1].split(';')[0]
                alt = info.split(f'{archaic_genome}_ALT=')[1].split(';')[0]
                gt = info.split(f'{archaic_genome}_GT=')[1].split(';')[0]
                if ancestral in [ref, alt] and '/' not in gt:
                    fields[f'{archaic_genome}_REF'].append(ref)
                    fields[f'{archaic_genome}_ALT'].append(alt)
                    fields[f'{archaic_genome}_GT'].append(int(gt))
                else:
                    fields[f'{archaic_genome}_REF'].append('-')
                    fields[f'{archaic_genome}_ALT'].append('-')
                    fields[f'{archaic_genome}_GT'].append(np.nan)
            except IndexError:
                fields[f'{archaic_genome}_REF'].append('-')
                fields[f'{archaic_genome}_ALT'].append('-')
                fields[f'{archaic_genome}_GT'].append(np.nan)

        fields['gnomad_ID'].append(get_basic_field_annotation(info, 'gnomad_ID', '-'))
        fields['segdup'].append(get_basic_field_annotation(info, 'gnomad_segdup', np.nan, return_int=True))
        fields['lcr'].append(get_basic_field_annotation(info, 'gnomad_lcr', np.nan, return_int=True))

        try:
            vep = info.split('gnomad_ensembl_vep=')[1].split(';')[0].split(',')
            csq = []
            for v in vep:
                if ancestral and ancestral != v.split('|')[0] and v.split('|')[0] in [ref, alt]:
                    csq.append(v.split('|')[1])
            if len(csq) == 0:
                fields['vep_consequence'].append(['-'])
            else:
                fields['vep_consequence'].append(csq)
        except IndexError:
            fields['vep_consequence'].append(['-'])

        fields['clinvar_ID'].append(get_basic_field_annotation(info, 'clinvar_ID', '-'))

        try:
            diseases = [dn for dn in info.split('clinvar_dn=')[1].split(';')[0].split('|')
                        if dn != 'not_provided' and dn != 'not_specified']
            if len(diseases) > 0:
                fields['clinvar_dn'].append(','.join(diseases))
            else:
                fields['clinvar_dn'].append('-')
        except IndexError:
            fields['clinvar_dn'].append('-')

        fields['clinvar_disdb'].append(get_basic_field_annotation(info, 'clinvar_disdb', '-'))
        fields['clinvar_sig'].append(get_basic_field_annotation(info, 'clinvar_sig', ['-'], return_list=True))
        fields['clinvar_mc'].append(get_basic_field_annotation(info, 'clinvar_mc', ['-'], return_list=True))

    # add annotations to dataframe
    for key, vals in fields.items():
        df[key] = vals
    return df


def parse_allele_frequencies(paths, pop):
    df = pd.concat([pd.read_csv(f, sep='\t', header=0) for f in paths])
    actual_ref = []
    position = []
    for id in df.ID.values:
        position.append(int(id.split(':')[1].split(',')[0][:-1]))
        actual_ref.append(id.split(':')[1].split(',')[0][-1])
    df['position'] = position
    df['actual_ref'] = actual_ref
    df.ALT_FREQS = np.where(df.ALT.values == df.actual_ref.values, 1 - df.ALT_FREQS.values, df.ALT_FREQS.values)
    df = df.loc[:, ['#CHROM', 'position', 'ALT_FREQS']]
    df.sort_values(['#CHROM', "position"], inplace=True)
    df.rename({"ALT_FREQS": f"AF_{pop}"}, inplace=True, axis=1)
    df['#CHROM'] = [f'chr{chrom}' for chrom in df['#CHROM'].values]
    df.set_index(['#CHROM', "position"], inplace=True, drop=True)
    return df


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Path to BED file annotated with basic stats and '
                                              'modern variant annotations. Modern variant annotations will be '
                                              'aggregated per putatively selected introgressed segment')
    parser.add_argument('-o', '--output', help='Path to write annotated BED file to. Modern variant annotations will '
                                               'be aggregated by putatively selected introgressed segment')
    parser.add_argument('-b', '--basic_stats', help='Path to BED file with putatively selected introgressed segments '
                                                    'annotated with basic stats. Used to infer column names')
    parser.add_argument('--header_modern_stats', help='Path to file containing header, i.e., column names, of '
                                                      'annotation fields for modern variants')
    parser.add_argument('--aa_af', nargs='+',
                        help='Path to files containing allele frequencies in African-American individuals')
    parser.add_argument('--afr_af', nargs='+', help='Path to files containing allele frequencies in African individuals')
    parser.add_argument('--eur_af', nargs='+',
                        help='Path to files containing allele frequencies in European individuals')
    parser.add_argument('--eas_af', nargs='+',
                        help='Path to files containing allele frequencies in East Asian individuals')
    parser.add_argument('--archaic_genomes', nargs='+', help='Archaic genomes that were used in the analysis')

    args = parser.parse_args()
    archaic_genomes = args.archaic_genomes

    # by default column operations will be mean for previously annotations
    column_operations = {}
    column_operations['AF_AA'] = concat_strings_in_aggregation_all
    column_operations['AF_AFR'] = concat_strings_in_aggregation_all
    column_operations['AF_EUR'] = concat_strings_in_aggregation_all
    column_operations['AF_EAS'] = concat_strings_in_aggregation_all
    # select columns to keep
    columns_to_keep = ['chrom', 'start', 'end', 'REF', 'ALT', "INFO", 'AF_AA', 'AF_AFR', 'AF_EUR', 'AF_EAS']
    # get additional column names from VCF
    columns = ['chrom', 'start', 'end']
    with open(args.header_modern_stats, 'r') as header:
        cols = header.readline().strip().split('\t')
    header.close()
    columns.extend(cols)
    columns.append('overlap')
    aa_allele_frequencies = parse_allele_frequencies(args.aa_af, 'AA')
    afr_allele_frequencies = parse_allele_frequencies(args.afr_af, 'AFR')
    eur_allele_frequencies = parse_allele_frequencies(args.eur_af, 'EUR')
    eas_allele_frequencies = parse_allele_frequencies(args.eas_af, 'EAS')
    # read data
    df = pd.read_csv(args.input, header=None, sep='\t', names=columns, na_values='.')
    # remove sites with missing annotations
    df = df[df['#CHROM'] != '.']
    df = df[~pd.isna(df['#CHROM'])].reset_index(drop=True)
    # join allele frequencies
    df = df.join(aa_allele_frequencies, on=['#CHROM', "POS"])
    df = df.join(afr_allele_frequencies, on=['#CHROM', "POS"])
    df = df.join(eur_allele_frequencies, on=['#CHROM', "POS"])
    df = df.join(eas_allele_frequencies, on=['#CHROM', "POS"])

    # keep only needed columns and throw out sites that are monomorphic in African-Americans
    df = df.loc[df.AF_AA > 0, columns_to_keep]
    df = parse_modern_variant_annotations(df, archaic_genomes)

    # calculate derived allele frequencies
    df.AF_AA = np.where(df.REF == df.ancestral, df.AF_AA, 1 - df.AF_AA)
    df.AF_AFR = np.where(df.REF == df.ancestral, df.AF_AFR, 1 - df.AF_AFR)
    df.AF_EUR = np.where(df.REF == df.ancestral, df.AF_EUR, 1 - df.AF_EUR)
    df.AF_EAS = np.where(df.REF == df.ancestral, df.AF_EAS, 1 - df.AF_EAS)

    # polarize archaic genotypes by ancestral allele state
    for archaic_genome in archaic_genomes:
        df[f'{archaic_genome}_GT'] = np.where(df[f'{archaic_genome}_REF'] != df.ancestral,
                                              2 - df[f'{archaic_genome}_GT'], df[f'{archaic_genome}_GT'])
        df[f'AFR_fixed_ancestral_{archaic_genome}_GT'] = np.where((df[f'{archaic_genome}_GT'] != 0) &
                                                                  (df.AF_AFR > 0),
                                                                  '',
                                                                  df[f'{archaic_genome}_GT'].astype(str))
        column_operations[f'{archaic_genome}_GT'] = concat_strings_in_aggregation_all
        column_operations[f'AFR_fixed_ancestral_{archaic_genome}_GT'] = concat_strings_in_aggregation_all
    # encode VEP annotations
    categories, encoding = onehotencoding_variant_annotations(df.vep_consequence.values)
    df[[f'vep_{cat}' for cat in categories]] = encoding
    for cat in categories:
        column_operations[f'vep_{cat}'] = np.sum

    # encode CLINVAR significances
    categories, encoding = onehotencoding_variant_annotations(df.clinvar_sig.values)
    df[[f'clinvar_sig_{cat}' for cat in categories]] = encoding
    for cat in categories:
        column_operations[f'clinvar_sig_{cat}'] = np.sum

    # encode CLINVAR molecular consequences
    categories, encoding = onehotencoding_variant_annotations(df.clinvar_mc.values)
    df[[f'clinvar_mc_{cat}' for cat in categories]] = encoding
    for cat in categories:
        column_operations[f'clinvar_mc_{cat}'] = np.sum

    # set aggregation operation
    column_operations['gnomad_ID'] = concat_strings_in_aggregation
    column_operations['clinvar_ID'] = concat_strings_in_aggregation
    column_operations['clinvar_dn'] = concat_strings_in_aggregation
    column_operations['clinvar_disdb'] = concat_strings_in_aggregation
    df.AF_AFR = df.AF_AFR.astype(str)
    df.AF_EUR = df.AF_EUR.astype(str)
    df.AF_EAS = df.AF_EAS.astype(str)
    df.AF_AA = df.AF_AA.astype(str)
    df.Altai_GT = df.Altai_GT.fillna('').astype(str)
    df['Vindija33.19_GT'] = df['Vindija33.19_GT'].fillna('').astype(str)
    df.Chagyrskaya_GT = df.Chagyrskaya_GT.fillna('').astype(str)
    df.Denisova_GT = df.Altai_GT.fillna('').astype(str)

    # aggregate
    df = df.groupby(['chrom', 'start', 'end']).agg(column_operations)

    basic_stats = pd.read_csv(args.basic_stats, header=0, sep='\t')
    # avoid nans
    columns_to_handle = ['B_distinct', 'gerp_distinct', 'phylop_distinct', 'phastcons_distinct', 'recomb_distinct',
                         'boosting_complete_eur_distinct', 'boosting_complete_eur_recent_distinct',
                         'boosting_complete_eur_ancient_distinct', 'boosting_complete_afr_distinct',
                         'boosting_complete_afr_recent_distinct', 'boosting_complete_afr_ancient_distinct',
                         'boosting_complete_eas_distinct', 'boosting_complete_eas_recent_distinct',
                         'boosting_complete_eas_ancient_distinct', 'boosting_incomplete_eur_distinct',
                         'boosting_incomplete_eur_recent_distinct', 'boosting_incomplete_eur_ancient_distinct',
                         'boosting_incomplete_afr_distinct', 'boosting_incomplete_afr_recent_distinct',
                         'boosting_incomplete_afr_ancient_distinct', 'boosting_incomplete_eas_distinct',
                         'boosting_incomplete_eas_recent_distinct', 'boosting_incomplete_eas_ancient_distinct']
    for col in columns_to_handle:
        basic_stats[col] = ['' if (isinstance(val, float) and np.isnan(val)) or val == '.'
                            else str(val) for val in basic_stats[col]]

    # join dataframes
    df = df.join(basic_stats.set_index(['chrom', 'start', 'end']))
    df.reset_index(inplace=True)
    # write to file
    df.chrom = [int(chrom.replace('chr', '')) for chrom in df.chrom]
    df.chrom = df.chrom.astype(int)
    df.start = df.start.astype(int)
    df.end = df.end.astype(int)
    df.to_csv(args.output, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
