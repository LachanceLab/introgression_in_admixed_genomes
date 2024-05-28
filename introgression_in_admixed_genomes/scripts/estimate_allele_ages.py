#!/usr/bin/env python3
import sys
import argparse
import tsinfer
import tskit
from cyvcf2 import VCF
import pandas as pd
import msprime
import tsdate
import numpy as np
import multiprocessing
import os
import allel


def get_sample_ids(file):
    """
    Reads file with sample IDs corresponding to a given population. One ID per line
    :param file: str, path to filename
    :return: list, sample IDs
    """
    samples = []
    with open(file, 'r') as iids:
        for line in iids:
            # exclude HGDP samples
            if not line.strip().startswith('HGDP'):
                samples.append(line.strip())
    return samples


def chromosome_length(vcf, chrom):
    """
    Get chromosome length of a provided cyvcf2 VCF object
    :param vcf: cyvc2 VCF object
    :param chrom: int, current chormosome numher
    :return: int, chromosome length of corresponding VCF
    """
    return vcf.seqlens[chrom - 1]


def get_archaic_variant(vcf, chrom, pos):
    """
    Return variants at specified position
    :param vcf: cyvcf2 VCF
    :param chrom: int, chromosome number
    :param pos: int, position of interest
    :return:
    """
    try:
        return [variant for variant in vcf(f'chr{chrom}:{pos}-{pos}')][0]
    except IndexError:
        return None


def get_alleles(alleles, variant):
    """
    Add alleles of focal variant to list of alleles if not already present
    :param alleles: list, list of alleles observed at the current position so far
    :param variant: cyvcf Variant, focal variant
    :return: list, updated list of unique alleles
    """
    if variant and len(variant.ALT) > 0:
        alleles.extend([v.upper() for v in variant.ALT if not v.upper() in alleles])
    return alleles


def get_genotypes(genotypes, variant, alleles):
    """
    Add genotypes of a focal variant to existing list of genotypes
    :param genotypes: array-like, list of genotypes per sample
    :param variant: cyvcf2 Variant, focal variant
    :param alleles: list, list of alleles that are observed at current position. Used to encode genotypes of focal variant
    :return: list, updated list of genotypes
    """
    if not variant or len(variant.ALT) == 0:
        genotypes.extend([0, 0])
    else:
        genotypes.extend([alleles.index(([variant.REF.upper()] +
                                         [v.upper() for v in variant.ALT])[g][0].upper())
                          for g in variant.genotypes[0][:2]])
    return genotypes


def split_vcf_into_chunks(vcf, threads):
    """
    Split variable sites in VCF into equal sized chunks
    :param vcf: str, path to input VCF
    :param threads: int, number CPUs
    :return: list, regions which split the chromosome into equal sized chunks in terms of number of variable sites
    """
    variable_sites = allel.read_vcf(vcf, fields=["variants/CHROM", "variants/POS"])
    positions = variable_sites["variants/POS"]
    chromosomes = variable_sites["variants/CHROM"]
    assert np.all(chromosomes == chromosomes[0])
    chrom = str(chromosomes[0])
    n_variants = len(positions)
    chunk_size = int(np.ceil(n_variants / threads))
    regions = []
    for part, i in enumerate(range(0, n_variants, chunk_size)):
        if part == threads - 1:
            regions.append(f'chr{chrom}:{positions[i]}-{positions[-1]}')
        else:
            regions.append(f'chr{chrom}:{positions[i]}-{positions[i + chunk_size - 1]}')
    return regions


def merge_sampledata_files(input_sampledata, output):
    """
    Merge multiple tsinfer.SampleData objects. This script assumes that all input_sampledata contain the exact same
    individuals/populations and have no overlapping sites. Furthermore, it assumes that the sites are only increasing
    with each input sample data.
    """
    # load sample files
    sample_files = []
    for cur_sample in input_sampledata:
        sample_files.append(tsinfer.load(cur_sample))
    pop_lookup = {}
    # intialize merged SampleData object
    with tsinfer.SampleData(sequence_length=sample_files[0].sequence_length) as merged:
        # add unique populations
        for i, population in enumerate(sample_files[0].populations()):
            pop_lookup[i] = merged.add_population(metadata=population.metadata)
        # add unique individuals
        for individual in sample_files[0].individuals():
            merged.add_individual(location=individual.location, metadata=individual.metadata,
                                  time=individual.time, flags=individual.flags,
                                  population=pop_lookup[individual.population],
                                  ploidy=len(individual.samples))
        # add variants from all input_sampledata objects in sequential order
        for i, sample in enumerate(sample_files):
            for variant in sample.variants():
                merged.add_site(position=variant.site.position,
                                genotypes=variant.genotypes[:sample_files[0].num_samples],
                                alleles=variant.alleles, metadata=variant.site.metadata,
                                time=variant.site.time, ancestral_allele=variant.site.ancestral_allele)
    # write to file and finalize
    merged_samples = merged.copy(output)
    merged_samples.finalise()
    return merged_samples


def add_diploid_sites(vcf_modern, vcf_altai, vcf_vindija, vcf_chagyrskaya, vcf_denisovan,
                      samples, region):
    """
    Read the sites in the vcf and add them to the samples object. Only add variants from archaic genomes for
    which modern humans are polymorphic.
    :param vcf_modern: cyvcf2 VCF, corresponding to modern individuals
    :param vcf_altai: cyvcf2 VCF, corresponding to Altai Neanderthal individual
    :param vcf_vindija: cyvcf2 VCF, corresponding to Vindija Neanderthal individual
    :param vcf_chagyrskaya: cyvcf2 VCF, corresponding to Chagyrskaya Neanderthal individual
    :param vcf_denisovan: cyvcf2 VCF, corresponding to Denisovan individual
    :param samples: tsinfer.SampleData, sample data to which to add sites
    :param region: str, region (chrom:start-end) for which to add sites from VCFs
    """
    pos = int(region.split(':')[1].split('-')[0]) - 1
    start = int(region.split(':')[1].split('-')[0])
    end = int(region.split(':')[1].split('-')[1])
    duplicate_positions = []

    for variant_modern in vcf_modern(region):  # Loop over variants, each assumed at a unique site
        if pos == variant_modern.POS:
            duplicate_positions.append(pos)
            continue
        elif variant_modern.POS < start:
            print(f"Site {variant_modern.POS} out of range")
            continue
        elif variant_modern.POS > end:
            print(f"Site {variant_modern.POS} out of range")
            break
        else:
            pos = variant_modern.POS
        alleles = [variant_modern.REF.upper()] + [v.upper() for v in variant_modern.ALT]
        # consider only biallelic sites
        if len(alleles) > 2:
            continue

        ancestral = variant_modern.INFO.get("ancestral", ".").upper()  # "." means unknown
        # skip sites for which we don't know the ancestral allele status
        if ancestral == "." or ancestral == "" or len(ancestral) > 1:
            continue
        # get archaic variants and ensure that reference allele match
        if vcf_altai:
            altai = get_archaic_variant(vcf_altai, variant_modern.CHROM, variant_modern.POS)
            if altai and altai.REF != variant_modern.REF:
                continue
        if vcf_vindija:
            vindija = get_archaic_variant(vcf_vindija, variant_modern.CHROM, variant_modern.POS)
            if vindija and vindija.REF != variant_modern.REF:
                continue
        if vcf_chagyrskaya:
            chagyrskaya = get_archaic_variant(vcf_chagyrskaya, variant_modern.CHROM, variant_modern.POS)
            if chagyrskaya and chagyrskaya.REF != variant_modern.REF:
                continue
        if vcf_denisovan:
            denisovan = get_archaic_variant(vcf_denisovan, variant_modern.CHROM, variant_modern.POS)
            if denisovan and denisovan.REF != variant_modern.REF:
                continue

        # get archaic alleles
        if vcf_altai:
            alleles = get_alleles(alleles, altai)
        if vcf_vindija:
            alleles = get_alleles(alleles, vindija)
        if vcf_chagyrskaya:
            alleles = get_alleles(alleles, chagyrskaya)
        if vcf_denisovan:
            alleles = get_alleles(alleles, denisovan)

        if ancestral not in alleles:
            alleles.append(ancestral)
        ancestral_allele = alleles.index(ancestral)

        # get genotypes
        genotypes = [g for row in variant_modern.genotypes for g in row[0:2]]

        if vcf_vindija:
            genotypes = get_genotypes(genotypes, vindija, alleles)
        if vcf_denisovan:
            genotypes = get_genotypes(genotypes, denisovan, alleles)
        if vcf_chagyrskaya:
            genotypes = get_genotypes(genotypes, chagyrskaya, alleles)
        if vcf_altai:
            genotypes = get_genotypes(genotypes, altai, alleles)
        samples.add_site(pos, genotypes, alleles, ancestral_allele=ancestral_allele)
    return duplicate_positions


def read_vcf_part(params):
    """
    Helper function to parse VCFs into tsinfer.SampleData objects in parallel
    :param params: tuple, list of arguments
    :return: str, path to file with partial SampleData
    """
    vcf_modern_path, vcf_altai_path, vcf_vindija_path, vcf_chagyrskaya_path, vcf_denisovan_path, \
    chrom, AFR_samples, EUR_samples, EAS_samples, all_samples, \
    generation_time, part, output_dir, region = params
    vcf_modern = VCF(vcf_modern_path, samples=all_samples)
    vcf_modern.set_index(f'{vcf_modern_path}.tbi')
    if vcf_altai_path:
        vcf_altai = VCF(vcf_altai_path)
        vcf_altai.set_index(f'{vcf_altai_path}.tbi')
    else:
        vcf_altai = None
    if vcf_vindija_path:
        vcf_vindija = VCF(vcf_vindija_path)
        vcf_vindija.set_index(f'{vcf_vindija_path}.tbi')
    else:
        vcf_vindija = None
    if vcf_chagyrskaya_path:
        vcf_chagyrskaya = VCF(vcf_chagyrskaya_path)
        vcf_chagyrskaya.set_index(f'{vcf_chagyrskaya_path}.tbi')
    else:
        vcf_chagyrskaya = None
    if vcf_denisovan_path:
        vcf_denisovan = VCF(vcf_denisovan_path)
        vcf_denisovan.set_index(f'{vcf_denisovan_path}.tbi')
    else:
        vcf_denisovan = None

    sequence_length = chromosome_length(vcf_modern, chrom)

    path = f"{output_dir}/partial_file_part{part}.samples"
    pop_lookup = {}
    with tsinfer.SampleData(path=f'{path}.tmp', sequence_length=sequence_length) as samples:
        # add populations
        pop_lookup['AFR'] = samples.add_population(metadata={'name': 'AFR_1KGP'})
        pop_lookup['EUR'] = samples.add_population(metadata={'name': 'EUR_1KGP'})
        pop_lookup['EAS'] = samples.add_population(metadata={'name': 'EAS_1KGP'})
        if vcf_altai_path:
            pop_lookup['altai'] = samples.add_population(metadata={'name': 'Nea_Altai'})
        if vcf_vindija_path:
            pop_lookup['vindija'] = samples.add_population(metadata={'name': 'Nea_Vindija33.19'})
        if vcf_chagyrskaya_path:
            pop_lookup['chagyrskaya'] = samples.add_population(metadata={'name': 'Nea_Chagyrskaya'})
        if vcf_denisovan_path:
            pop_lookup['denisovan'] = samples.add_population(metadata={'name': 'Denisova'})
        # add individuals
        for sample in [sample for sample in vcf_modern.samples if sample in AFR_samples]:
            samples.add_individual(ploidy=2, metadata={'iid': sample},
                                   population=pop_lookup['AFR'], time=0)
        for sample in [sample for sample in vcf_modern.samples if sample in EUR_samples]:
            samples.add_individual(ploidy=2, metadata={'iid': sample},
                                   population=pop_lookup['EUR'], time=0)
        for sample in [sample for sample in vcf_modern.samples if sample in EAS_samples]:
            samples.add_individual(ploidy=2, metadata={'iid': sample},
                                   population=pop_lookup['EAS'], time=0)
        if vcf_vindija_path:
            samples.add_individual(ploidy=2, metadata={'iid': 'vindija'}, population=pop_lookup['vindija'],
                                   time=50000 / generation_time)
        if vcf_denisovan_path:
            samples.add_individual(ploidy=2, metadata={'iid': 'denisovan'}, population=pop_lookup['denisovan'],
                                   time=63900 / generation_time)
        if vcf_chagyrskaya_path:
            samples.add_individual(ploidy=2, metadata={'iid': 'chagyrskaya'},
                                   population=pop_lookup['chagyrskaya'], time=80000 / generation_time)
        if vcf_altai_path:
            samples.add_individual(ploidy=2, metadata={'iid': 'altai'}, population=pop_lookup['altai'],
                                   time=110000 / generation_time)
        # add sites
        duplicate_positions = add_diploid_sites(vcf_modern, vcf_altai, vcf_vindija, vcf_chagyrskaya, vcf_denisovan,
                                                samples, region)
    positions = samples.sites_position[:]
    sites = np.array([s.id for s in samples.sites()])
    to_keep = sites[~np.isin(positions, duplicate_positions)]
    filtered_sites = samples.subset(sites=to_keep)
    filtered_sites.copy(path=path)
    filtered_sites.close()
    samples.close()
    os.remove(f'{path}.tmp')
    return path


def get_chromosome_arm_coords(args):
    chromosome_arms = pd.read_csv(args.chromosome_arms_df, sep='\t', names=['chrom', 'start', 'end', 'arm'])
    chromosome_arms = chromosome_arms[(chromosome_arms.chrom == args.chromosome) &
                                      (chromosome_arms.arm == args.chromosome_arm)]
    start = chromosome_arms.start.values[0]
    end = chromosome_arms.end.values[0]
    return start, end


def get_recombination_map(hapmap_path, samples, start, end):
    rate_map = msprime.RateMap.read_hapmap(hapmap_path, sequence_length=samples.sequence_length)
    # substitute NaNs and 0s
    rates = np.nan_to_num(rate_map.rate, nan=1e-20)
    rates[rates == 0] = 1e-20
    rate_map = msprime.RateMap(position=np.concatenate([rate_map.left, [rate_map.right[-1]]]), rate=rates)
    rate_map = rate_map.slice(start, min([end, samples.sequence_length]))
    return rate_map


def format_outputdir(output_dir):
    if output_dir and not os.path.isdir(output_dir):
        os.makedirs(output_dir)
    if output_dir and not output_dir.endswith('/'):
        output_dir += '/'
    return output_dir


def parse_vcf(args):
    output_dir = format_outputdir(args.output_dir)
    AFR_samples = get_sample_ids(args.afr_ids)
    EUR_samples = get_sample_ids(args.eur_ids)
    EAS_samples = get_sample_ids(args.eas_ids)
    all_samples = AFR_samples.copy()
    all_samples.extend(EUR_samples)
    all_samples.extend(EAS_samples)
    regions = split_vcf_into_chunks(args.vcf_modern, args.threads)
    ready_to_map = [(args.vcf_modern, args.vcf_altai, args.vcf_vindija,
                     args.vcf_chagyrskaya, args.vcf_denisovan, args.chromosome,
                     AFR_samples, EUR_samples, EAS_samples,
                     all_samples, args.generation_time, i, output_dir, region)
                    for i, region in enumerate(regions)]
    pool = multiprocessing.Pool(processes=args.threads)
    pool.map(read_vcf_part, ready_to_map)
    pool.close()
    pool.join()


def merge_sampledata(args):
    output_dir = format_outputdir(args.output_dir)
    files_part_mapping = {int(fn.split('_part')[1].split('.')[0]): fn for fn in args.sampledata}
    sample_files = [files_part_mapping[i] for i in range(len(args.sampledata))]
    merge_sampledata_files(sample_files, f'{output_dir}chr{args.chromosome}.samples')


def infer_chromosome_arm(args):
    output_dir = format_outputdir(args.output_dir)
    # get start end coordinates of chromosome arm
    start, end = get_chromosome_arm_coords(args)
    # load samples
    samples = tsinfer.load(args.sampledata)
    positions = samples.sites_position[:]
    sites = np.array([s.id for s in samples.sites()])
    # tsinfer takes site positions when excluding them
    sites_to_exclude = positions[np.where((positions < start) | (positions >= end))[0]]
    # subset takes sites IDs
    sites_to_keep = sites[np.where((positions >= start) & (positions < end))[0]]
    # subset sites to chromosome arm
    arm = samples.subset(sites=sites_to_keep)
    arm_samples = arm.copy(path=args.sampledata.replace('.samples', f'_{args.chromosome_arm}_arm.samples'))
    arm_samples.finalise()
    # subset to modern individuals
    modern_samples = arm.subset(np.where(arm.individuals_time[:] == 0)[0])
    # load ratemap
    rate_map = get_recombination_map(args.hapmap, samples, start, end)
    # infer tree sequence
    inferred_ts = tsinfer.infer(modern_samples, recombination_rate=rate_map.mean_rate,
                                num_threads=args.threads, exclude_positions=sites_to_exclude)
    inferred_ts.dump(f'{output_dir}chr{args.chromosome}_{args.chromosome_arm}_arm.trees')


def date_alleles(args):
    ts = tskit.load(args.trees)
    # Removes unary nodes (currently required in tsdate), keeps historical-only sites
    preprocessed_ts = tsdate.preprocess_ts(ts, filter_sites=False)
    dated_ts = tsdate.date(preprocessed_ts, Ne=args.ne, mutation_rate=args.mutation_rate,
                           num_threads=args.threads)
    dated_ts.dump(args.trees.replace('.trees', '_dated.trees'))


def get_dated_sampledata_from_ts(args):
    # load trees
    ts = tskit.load(args.trees)
    # load samples
    samples = tsinfer.load(args.sampledata)
    # get times
    sites_time = tsdate.sites_time_from_ts(ts)
    # update sample times
    dated_samples = tsdate.add_sampledata_times(samples, sites_time)
    # save
    dated_s = dated_samples.copy(args.trees.replace('.trees', '.samples'))
    dated_s.finalise()
    dated_s.close()


def generate_ancestors_ts(args):
    output_dir = format_outputdir(args.output_dir)
    # load samples
    samples = tsinfer.load(args.sampledata)
    start, end = get_chromosome_arm_coords(args)
    # get positions to exclude
    positions = samples.sites_position[:]
    sites_to_exclude = positions[np.where((positions < start) | (positions >= end))[0]]
    # load ratemap
    rate_map = get_recombination_map(args.hapmap, samples, start, end)
    ancestors = tsinfer.generate_ancestors(samples, exclude_positions=sites_to_exclude)
    ancestors_w_proxy = ancestors.insert_proxy_samples(samples, allow_mutation=True)
    ancestors_ts = tsinfer.match_ancestors(samples, ancestors_w_proxy,
                                           recombination_rate=rate_map.mean_rate,
                                           num_threads=args.threads, path_compression=False)
    ancestors_ts.dump(f'{output_dir}chr{args.chromosome}_{args.chromosome_arm}_arm_ancestors.trees')


def reinfer_ts(args):
    output_dir = format_outputdir(args.output_dir)
    ts = tskit.load(args.trees)
    samples = tsinfer.load(args.sampledata)
    # get chromosome arm coords
    start, end = get_chromosome_arm_coords(args)
    # load ratemap
    rate_map = get_recombination_map(args.hapmap, samples, start, end)
    reinferred_ts = tsinfer.match_samples(samples, ts, force_sample_times=True,
                                          recombination_rate=rate_map.mean_rate, num_threads=args.threads)
    reinferred_ts.dump(f"{output_dir}chr{args.chromosome}_{args.chromosome_arm}_arm_reinferred.trees")


def get_allele_ages(args):
    ts = tskit.load(args.trees)
    samples = tsinfer.load(args.sampledata)

    populations = samples.individuals_population[samples.samples_individual[:]]
    ancestral_alleles = samples.sites_ancestral_allele
    # polarize genotypes by ancestral allele status
    genotypes = samples.sites_genotypes[:]
    polarized_genotypes = np.zeros_like(genotypes)
    for i in range(genotypes.shape[0]):
        polarized_genotypes[i, :] = np.where(genotypes[i] == ancestral_alleles[i], 0, 1)
    # get TMRCAS for each site in kya
    sites_tmrca = tsdate.sites_time_from_ts(ts, unconstrained=False, node_selection='arithmetic') *\
                  args.generation_time / 1000
    # get current positions
    positions = samples.sites_position[:]
    df_dict = {'CHROM': np.repeat(args.chromosome, positions.shape[0]),
               'POS': positions.astype(int),
               'TMRCA_kya': sites_tmrca}
    # calculate population-wise allele frequencies
    for pop in samples.populations():
        df_dict[pop.metadata['name']] = polarized_genotypes[:, populations == pop.id].mean(axis=1)
    df = pd.DataFrame.from_dict(df_dict)
    neanderthals = [col for col in df.columns if col.startswith('Nea_')]
    df.All_Nea = df.loc[:, neanderthals].mean(axis=1).values
    df.to_csv(f'{args.trees.replace("reinferred.trees", "tmrca_summarized.tab")}', sep='\t', header=True,
              index=False)


def main(argv):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()

    subparser = subparsers.add_parser('parse_vcfs')
    subparser.add_argument('--vcf_modern', help='Modern VCF filename')
    subparser.add_argument('--vcf_altai', help='Altai VCF filenames', required=False, default=None)
    subparser.add_argument('--vcf_vindija', help='Vindija VCF filenames', required=False, default=None)
    subparser.add_argument('--vcf_chagyrskaya', help='Chagyrskaya VCF filenames', required=False, default=None)
    subparser.add_argument('--vcf_denisovan', help='Denisovan VCF filename', required=False, default=None)
    subparser.add_argument('--afr_ids', help='File containing sample IDs of African individuals')
    subparser.add_argument('--eur_ids', help='File containing sample IDs of European individuals')
    subparser.add_argument('--eas_ids', help='File containing sample IDs of East-Asian individuals')
    subparser.add_argument('-o', '--output_dir', help='Output directory')
    subparser.add_argument('--chromosome', type=int, help='Current chromosome number')
    subparser.add_argument('-g', '--generation_time', type=int, help='Assumed generation time for dating in years [29]',
                           default=29)
    subparser.add_argument('-t', '--threads', type=int, help='Number of CPUs to use')
    subparser.set_defaults(func=parse_vcf)

    subparser = subparsers.add_parser('merge_sampledata')
    subparser.add_argument('--sampledata', nargs='+', help='tsinfer.SampleData file(s). '
                                                           'You can specify multiple when merging')
    subparser.add_argument('-o', '--output_dir', help='Output directory')
    subparser.add_argument('--chromosome', type=int, help='Current chromosome number')
    subparser.set_defaults(func=merge_sampledata)

    subparser = subparsers.add_parser('infer_chromosome_arm')
    subparser.add_argument('--sampledata', help='tsinfer.SampleData file')
    subparser.add_argument('--chromosome_arms_df',
                           help='Filename of BED file containing coordinates of chromosome arms')
    subparser.add_argument('--hapmap', help='Filename of HapMap recombination names. '
                                            'Must be readable with msprime.RateHap.read_hapmap')
    subparser.add_argument('-o', '--output_dir', help='Output directory')
    subparser.add_argument('--chromosome_arm', help='Current chromosome arm to consider (q or p)')
    subparser.add_argument('--chromosome', type=int, help='Current chromosome number')
    subparser.add_argument('-t', '--threads', type=int, help='Number of CPUs to use')
    subparser.set_defaults(func=infer_chromosome_arm)


    subparser = subparsers.add_parser('date_alleles')
    subparser.add_argument('--trees', help='Input TreeSequence')
    subparser.add_argument('--ne', type=int, help='Effective population size passed to tsdate [10000]', default=10000)
    subparser.add_argument('-m', '--mutation_rate', type=float, help='Assumed mutation rate [1.25e-8]', default=1.25e-8)
    subparser.add_argument('-t', '--threads', type=int, help='Number of CPUs to use')
    subparser.set_defaults(func=date_alleles)

    subparser = subparsers.add_parser('get_dated_sampledata_from_ts')
    subparser.add_argument('--sampledata', help='tsinfer.SampleData file')
    subparser.add_argument('--trees', help='Input TreeSequence')
    subparser.set_defaults(func=get_dated_sampledata_from_ts)

    subparser = subparsers.add_parser('generate_ancestors_ts')
    subparser.add_argument('--sampledata', help='tsinfer.SampleData file(s).'
                                                           'You can specify multiple when merging')
    subparser.add_argument('--chromosome_arms_df',
                           help='Filename of BED file containing coordinates of chromosome arms')
    subparser.add_argument('--hapmap', help='Filename of HapMap recombination names. '
                                         'Must be readable with msprime.RateHap.read_hapmap')
    subparser.add_argument('-o', '--output_dir', help='Output directory')
    subparser.add_argument('--chromosome_arm', help='Current chromosome arm to consider (q or p)')
    subparser.add_argument('--chromosome', type=int, help='Current chromosome number')
    subparser.add_argument('-t', '--threads', type=int, help='Number of CPUs to use')
    subparser.set_defaults(func=generate_ancestors_ts)

    subparser = subparsers.add_parser('reinfer_ts')
    subparser.add_argument('--trees', help='Input TreeSequence')
    subparser.add_argument('--sampledata', help='tsinfer.SampleData file')
    subparser.add_argument('--chromosome_arms_df',
                           help='Filename of BED file containing coordinates of chromosome arms')
    subparser.add_argument('--hapmap', help='Filename of HapMap recombination names. '
                                            'Must be readable with msprime.RateHap.read_hapmap')
    subparser.add_argument('-o', '--output_dir', help='Output directory')
    subparser.add_argument('--chromosome_arm', help='Current chromosome arm to consider (q or p)')
    subparser.add_argument('--chromosome', type=int, help='Current chromosome number')
    subparser.add_argument('-t', '--threads', type=int, help='Number of CPUs to use')
    subparser.set_defaults(func=reinfer_ts)

    subparser = subparsers.add_parser('get_allele_ages')
    subparser.add_argument('--trees', help='Input TreeSequence')
    subparser.add_argument('--sampledata', help='tsinfer.SampleData file')
    subparser.add_argument('-g', '--generation_time', type=int, help='Assumed generation time for dating in years [29]',
                           default=29)
    subparser.add_argument('--chromosome', type=int, help='Current chromosome number')
    subparser.set_defaults(func=get_allele_ages)

    args = parser.parse_args()
    args.func(args)


if __name__ == '__main__':
    main(sys.argv[1:])
