#!/usr/bin/env python3
import sys
import argparse
import stdpopsim
import numpy as np
import msprime
import math
import subprocess


def run_simulations(samples, generation_time=25, mutation_rate=2.36e-8, simple=False):
    """
    Simulate simple human evolution with archaic admixture and an admixed American population
    :param samples: dict, target sample sizes
    :param generation_time: int
    :param mutation_rate: float
    :param simple, boolean, whether to run Gravel or simple demographic model
    :return: list, list: list of treesequences of simulated samples, list of masks
    """
    if simple:
        demography = msprime.Demography()
        demography.add_population(name='AFR', initial_size=20000, initially_active=True)
        demography.add_population(name='EUR', initial_size=10000)
        demography.add_population(name='EAS', initial_size=10000)
        demography.add_population(name='AMR', initial_size=30000)
        # ancestral Eurasian population
        demography.add_population(name='EUR_EAS', initial_size=10000)
        # Altai Neanderthal
        demography.add_population(name='NEA', initial_size=2800, default_sampling_time=110000 / generation_time)
        # Denisovan
        demography.add_population(name='DEN', initial_size=2600, default_sampling_time=63900 / generation_time)
        demography.add_population(
            name="AMH", description="Anatomically modern humans", initial_size=20000
        )
        demography.add_population(name='NEA_DEN', initial_size=2800)
        demography.add_population(name='ANC', initial_size=20000)

        # American admixture 15 generations ago
        demography.add_admixture(time=15, derived="AMR", ancestral=["AFR", "EUR", "EAS"],
                                 proportions=[0.8, 0.19, 0.01])
        # Eurasian split
        demography.add_population_split(time=920, derived=['EUR', "EAS"], ancestral="EUR_EAS")
        # Neanderthal introgression 50 kya
        demography.add_mass_migration(time=1500, source="EUR_EAS", dest="NEA", proportion=0.05)
        # OOA ~ 51 kya
        demography.add_population_split(time=2040, derived=["EUR_EAS", "AFR"], ancestral="AMH")
        # increase population size ~150 kya assuming 25 y/g
        demography.add_population_split(time=5920, derived=['AMH'], ancestral='ANC')
        # Neanderthal - Denisovan split ~ 500kya
        demography.add_population_split(time=20000, derived=["DEN", "NEA"], ancestral="NEA_DEN")
        # Neanderthal/Denisovan split from AMH ~700kya
        demography.add_population_split(time=28000, derived=["NEA_DEN"], ancestral="ANC")

    else:
        # Gravel's model with American admixture by Browning and
        demography = msprime.Demography()
        demography.add_population(name='AFR', initial_size=14474, initially_active=True)
        demography.add_population(name='EUR', initial_size=1032 * np.exp(3.8e-3 * 920), growth_rate=3.8e-3)
        demography.add_population(name='EAS', initial_size=554 * np.exp(4.8e-3 * 920), growth_rate=4.8e-3)
        demography.add_population(name='AMR', initial_size=30000 * np.exp(0.05 * 15), growth_rate=0.05)
        # ancestral Eurasian population
        demography.add_population(name='EUR_EAS', initial_size=1861)
        # Altai Neanderthal
        demography.add_population(name='NEA', initial_size=2800, default_sampling_time=110000 / generation_time)
        # Denisovan
        demography.add_population(name='DEN', initial_size=2600, default_sampling_time=63900 / generation_time)
        demography.add_population(
            name="AMH", description="Anatomically modern humans", initial_size=14474
        )
        demography.add_population(name='NEA_DEN', initial_size=2800)
        demography.add_population(name='ANC', initial_size=7310)

        # modern migration
        demography.set_symmetric_migration_rate(["EAS", "EUR"], rate=3.11e-5)
        demography.set_symmetric_migration_rate(["AFR", "EUR"], rate=2.5e-5)
        demography.set_symmetric_migration_rate(["AFR", "EAS"], rate=7.8e-6)

        # American admixture 15 generations ago
        demography.add_admixture(time=15, derived="AMR", ancestral=["AFR", "EUR", "EAS"],
                                 proportions=[0.8, 0.19, 0.01])
        # Eurasian split
        demography.add_population_split(time=920, derived=['EUR', "EAS"], ancestral="EUR_EAS")
        # migration between AFR and OOA
        demography.add_symmetric_migration_rate_change(time=920, populations=["AFR", "EUR_EAS"], rate=1.5e-4)

        # Neanderthal introgression 50 kya
        demography.add_mass_migration(time=1500, source="EUR_EAS", dest="NEA", proportion=0.02)
        # OOA ~ 51 kya
        demography.add_population_split(time=2040, derived=["EUR_EAS", "AFR"], ancestral="AMH")
        # increase population size ~150 kya assuming 25 y/g
        demography.add_population_split(time=5920, derived=['AMH'], ancestral='ANC')
        # Neanderthal - Denisovan split ~ 500kya
        demography.add_population_split(time=20000, derived=["DEN", "NEA"], ancestral="NEA_DEN")
        # Neanderthal/Denisovan split from AMH ~700kya
        demography.add_population_split(time=28000, derived=["NEA_DEN"], ancestral="ANC")

    # simulate 10 chromosomes using chr16 recombination map --> ~900Mb
    species = stdpopsim.get_species("HomSap")
    contigs = [species.get_contig(f'chr16', genetic_map='HapMapII_GRCh38', mutation_rate=mutation_rate)
               for chrom in np.arange(1, 11)]
    # this regions as very low recombination rates and leds to very long Neanderthal introgressed segments
    masks = {f'chr{i}': [31000000, 47000000] for i, contig in enumerate(contigs, 1)}
    # concatenate chromosomes
    rates = contigs[0].recombination_map.rate
    positions = np.concatenate([contigs[0].recombination_map.left,
                                [contigs[0].recombination_map.right[-1]]])
    for contig in contigs[1:]:
        positions = np.concatenate([positions, contig.recombination_map.left + positions[-1] + 1,
                                    [contig.recombination_map.right[-1] + 1 + positions[-1]]])
        rates = np.concatenate([rates, [math.log(2)], contig.recombination_map.rate])
    # get nan intervals
    rate_map = msprime.RateMap(position=positions, rate=rates)
    nan_intervals = [[left, right] for left, right in zip(rate_map.left[np.isnan(rate_map.rate)],
                                                          rate_map.right[np.isnan(rate_map.rate)])]
    # fill in nans
    rates_filled_nan = np.nan_to_num(rates, nan=0)
    rate_map = msprime.RateMap(position=positions, rate=rates_filled_nan)
    ts = msprime.sim_ancestry(samples=samples, demography=demography, recombination_rate=rate_map,
                              model='dtwf')
    # delete nan intervals
    ts = ts.delete_intervals(nan_intervals, simplify=False)
    # simulate mutations
    ts = msprime.sim_mutations(ts, rate=mutation_rate)
    # split chromosomes
    ts_chroms = []
    start = 0
    for contig in contigs:
        end = start + contig.recombination_map.right[-1]

        chrom_ts = ts.keep_intervals([[start, end]], simplify=False).trim()
        ts_chroms.append(chrom_ts)
        start += contig.recombination_map.right[-1] + 1
    return ts_chroms, masks


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--n_afr', type=int, help='Sample size for YRI (representative for African ancestry) [1067]',
                        default=1067)
    parser.add_argument('--n_afr_ref', type=int, default=504,
                        help='Sample size for YRI for PCA (representative for African ancestry) [504]')
    parser.add_argument('--n_eur', type=int, help='Sample size for CEU (representative for European ancestry) [10503]',
                        default=10503)
    parser.add_argument('--n_eur_ref', type=int, help='Sample size for CEU reference for PCA'
                                                      '(representative for European ancestry) [503]',
                        default=503)
    parser.add_argument('--n_eas', type=int, help='Sample size for CHB (representative for East Asian ancestry) [575]',
                        default=575)
    parser.add_argument('--n_eas_ref', type=int, default=504,
                        help='Sample size for CHB for PCA (representative for East Asian ancestry) [504]')
    parser.add_argument('--n_amr', type=int, help='Sample size for AMR (representative for African-American ancestry)')
    parser.add_argument('--trees', help='File to which to save the tree sequence. '
                                        'Must allow string formatting with chrom')
    parser.add_argument('--vcf', help='VCF file to which to write variants. Must allow string formatting with chrom')
    parser.add_argument('--ref', help='Reference panel. Will contain two columns: IID & FID. '
                                      'Must allow string formatting with chrom')
    parser.add_argument('--aa', help='Sample IDs of admixed individuals. Must allow string formatting with chrom')
    parser.add_argument('--afr', help='Sample IDs of African individuals. Must allow string formatting with chrom')
    parser.add_argument('--eur', help='Sample IDs of European individuals. Must allow string formatting with chrom')
    parser.add_argument('--eas', help='Sample IDs of East Asian individuals. Must allow string formatting with chrom')
    parser.add_argument('--nea', help='Sample ID of Neanderthal individual. Must allow string formatting with chrom')
    parser.add_argument('--den', help='Sample ID of Denisovan individual. Must allow string formatting with chrom')
    parser.add_argument('-g', '--genomefile', help='Genome file name. Must allow string formatting with chrom')
    parser.add_argument('--masks', help='Masked regions of chromosome. Must allow string formatting with chrom.')
    parser.add_argument('--simple', action='store_true', help='Whether to run full or simple demographic model.',
                        default=False)
    args = parser.parse_args()
    sample_sizes = [msprime.SampleSet(num_samples=args.n_afr, population='AFR'),
                    msprime.SampleSet(num_samples=args.n_eur, population='EUR'),
                    msprime.SampleSet(num_samples=args.n_eas, population='EAS'),
                    msprime.SampleSet(num_samples=min([63510, args.n_amr]), population='AMR'),
                    msprime.SampleSet(num_samples=1, population='NEA'),
                    msprime.SampleSet(num_samples=1, population='DEN')]
    ts_chroms, masks = run_simulations(sample_sizes, simple=args.simple)
    for chrom, ts in enumerate(ts_chroms, 1):
        ts.dump(args.trees.format(chrom=chrom))
        with open(args.masks.format(chrom=chrom), 'w') as mf:
            c_mask = masks[f'chr{chrom}']
            mf.write(f'{chrom}\t{c_mask[0]}\t{c_mask[1]}\n')
        mf.close()
    for chrom, ts in enumerate(ts_chroms, 1):
        with open(args.vcf.format(chrom=chrom), "w") as vcf_file:
            ts.write_vcf(vcf_file, contig_id=chrom)
        vcf_file.close()
        subprocess.run(["bgzip",  args.vcf.format(chrom=chrom)])
        # generate eur ref panel only once --> to ensure consistence across chromosomes
        if chrom == 1:
            eur_ref_panel = np.random.choice([indv.id for indv in ts.individuals()
                                              if ts.node(indv.nodes[0]).population == 1 and
                                              ts.node(indv.nodes[0]).time == 0],
                                             args.n_eur_ref, replace=False)
            afr_ref_panel = np.random.choice([indv.id for indv in ts.individuals()
                                              if ts.node(indv.nodes[0]).population == 0 and
                                              ts.node(indv.nodes[0]).time == 0], args.n_afr_ref, replace=False)
            eas_ref_panel = np.random.choice([indv.id for indv in ts.individuals()
                                              if ts.node(indv.nodes[0]).population == 2 and
                                              ts.node(indv.nodes[0]).time == 0], args.n_eas_ref, replace=False)
            aa_indvs = open(args.aa, 'w')
            afr_indvs = open(args.afr, 'w')
            eur_indvs = open(args.eur, 'w')
            eas_indvs = open(args.eas, 'w')
            nea_indv = open(args.nea, 'w')
            den_indv = open(args.den, 'w')
            ref_panel = open(args.ref, 'w')
            for indv in ts.individuals():
                if ts.node(indv.nodes[0]).population == 0:
                    if indv.id in afr_ref_panel:
                        ref_panel.write(f'tsk_{indv.id}\tAFR\n')
                    afr_indvs.write(f'tsk_{indv.id}\n')
                elif ts.node(indv.nodes[0]).population == 1:
                    if indv.id in eur_ref_panel:
                        ref_panel.write(f'tsk_{indv.id}\tEUR\n')
                    eur_indvs.write(f'tsk_{indv.id}\n')
                elif ts.node(indv.nodes[0]).population == 2:
                    if indv.id in eas_ref_panel:
                        ref_panel.write(f'tsk_{indv.id}\tEAS\n')
                    eas_indvs.write(f'tsk_{indv.id}\n')
                elif ts.node(indv.nodes[0]).population == 3:
                    aa_indvs.write(f'tsk_{indv.id}\n')
                elif ts.node(indv.nodes[0]).population == 5:
                    nea_indv.write(f'tsk_{indv.id}\n')
                elif ts.node(indv.nodes[0]).population == 6:
                    den_indv.write(f'tsk_{indv.id}\n')
            ref_panel.close()
            nea_indv.close()
            den_indv.close()
            afr_indvs.close()
            eas_indvs.close()
            eur_indvs.close()
            aa_indvs.close()
        with open(args.genomefile.format(chrom=chrom), 'w') as gf:
            gf.write(f"{chrom}\t{int(ts.sequence_length)}\n")
        gf.close()

    # # sample F1 generation
    # for chrom, ts in enumerate(ts_chroms, 1):
    #     # split tree sequence by populations
    #     f1_nodes = [node for node in ts.samples() if ts.node(node).time != 0]
    #     ts_f1 = ts.simplify(f1_nodes)
    #     with open(args.vcf.format(chrom=chrom).replace('.vcf', '_f1.vcf'), 'w') as vcf_file:
    #         ts_f1.write_vcf(vcf_file, contig_id=chrom)
    #     vcf_file.close()
    #     ref_panel = open(args.ref.format(chrom=chrom).replace('.txt', '_f1.txt'), 'w')
    #     # generate eur ref panel only once --> to ensure consistence across chromosomes
    #     if chrom == 1:
    #         eur_ref_panel = np.random.choice([indv.id for indv in ts_f1.individuals()
    #                                           if ts_f1.node(indv.nodes[0]).population == 1 and ts_f1.node(indv.nodes[0]).time == 14],
    #                                          args.n_eur_ref, replace=False)
    #     afr_indvs = open(args.afr.format(chrom=chrom).replace('.txt', '_f1.txt'), 'w')
    #     eur_indvs = open(args.eur.format(chrom=chrom).replace('.txt', '_f1.txt'), 'w')
    #     eas_indvs = open(args.eas.format(chrom=chrom).replace('.txt', '_f1.txt'), 'w')
    #     nea_indv = open(args.nea.format(chrom=chrom).replace('.txt', '_f1.txt'), 'w')
    #     den_indv = open(args.den.format(chrom=chrom).replace('.txt', '_f1.txt'), 'w')
    #     for indv in ts_f1.individuals():
    #         if ts_f1.node(indv.nodes[0]).population == 0:
    #             ref_panel.write(f'tsk_{indv.id}\tAFR\n')
    #             afr_indvs.write(f'tsk_{indv.id}\n')
    #         elif ts_f1.node(indv.nodes[0]).population == 1:
    #             if indv.id in eur_ref_panel:
    #                 ref_panel.write(f'tsk_{indv.id}\tEUR\n')
    #             eur_indvs.write(f'tsk_{indv.id}\n')
    #         elif ts_f1.node(indv.nodes[0]).population == 2:
    #             ref_panel.write(f'tsk_{indv.id}\tEAS\n')
    #             eas_indvs.write(f'tsk_{indv.id}\n')
    #         elif ts_f1.node(indv.nodes[0]).population == 5:
    #             nea_indv.write(f'tsk_{indv.id}\n')
    #         elif ts_f1.node(indv.nodes[0]).population == 6:
    #             den_indv.write(f'tsk_{indv.id}\n')
    #     ref_panel.close()
    #     nea_indv.close()
    #     den_indv.close()
    #     afr_indvs.close()
    #     eas_indvs.close()
    #     eur_indvs.close()
    #     with open(args.genomefile.format(chrom=chrom).replace('.bed', '_f1.bed'), 'w') as gf:
    #         gf.write(f"chr{chrom}\t{int(ts_f1.sequence_length)}\n")
    #     gf.close()


if __name__ == '__main__':
    main(sys.argv[1:])
