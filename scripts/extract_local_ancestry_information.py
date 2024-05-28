#!/usr/bin/env python3
import gzip
import sys
import argparse
import pandas as pd


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--lai', help='VCF with local ancestry information as outputted by FLARE')
    parser.add_argument('--output_phase0', help='BED file with local ancestry information for phase 0')
    parser.add_argument('--output_phase1', help='BED file with local ancestry information for phase 1')
    args = parser.parse_args()
    # open input file
    vcf_file = gzip.open(args.lai, 'rt')
    data_phase_0 = {}
    data_phase_1 = {}
    for line in vcf_file:
        if line.startswith('##') and not line.startswith('##ANCESTRY'):
            continue
        # parse ancestry encoding
        elif line.startswith("##ANCESTRY"):
            id_ancestry_mapping = {anc.split('=')[1]: anc.split('=')[0] for anc in
                                   line.strip().replace('##ANCESTRY=<', '').replace('>', '').split(',')}
        elif line.startswith('#CHROM'):
            # get sample names
            samples = line.strip().split('\t')[9:]
        else:
            line = line.strip().split('\t')
            chrom = line[0]
            pos = int(line[1])
            # iterate through lines
            for sample, lai in zip(samples, line[9:]):
                if not sample in data_phase_0:
                    # initialize dictionary for each phase
                    data_phase_0[sample] = dict()
                    data_phase_1[sample] = dict()
                    # start of current haplotype for each individual
                    data_phase_0[sample]['start'] = pos
                    data_phase_1[sample]['start'] = pos
                    # pointer to previous data point for each individual
                    data_phase_0[sample]['prev_pos'] = pos
                    data_phase_1[sample]['prev_pos'] = pos
                    # local ancestry of each individual
                    data_phase_0[sample]['c_ancestry'] = lai.split(':')[1]
                    data_phase_1[sample]['c_ancestry'] = lai.split(':')[2]
                    # initialize empty lists to store data
                    data_phase_0['chrom'] = []
                    data_phase_1['chrom'] = []
                    data_phase_0['start'] = []
                    data_phase_1['start'] = []
                    data_phase_0['end'] = []
                    data_phase_1['end'] = []
                    data_phase_0['la'] = []
                    data_phase_1['la'] = []
                    data_phase_0['IID'] = []
                    data_phase_1['IID'] = []
                    continue

                # switch in ancestry --> store haplotype info --> chrom, start, end, ancestry, IID
                if data_phase_0[sample]['c_ancestry'] != lai.split(':')[1]:
                    data_phase_0['chrom'].append(chrom)
                    data_phase_0['start'].append(data_phase_0[sample]['start'] - 1)
                    data_phase_0['end'].append(data_phase_0[sample]['prev_pos'])
                    data_phase_0['la'].append(id_ancestry_mapping[data_phase_0[sample]['c_ancestry']])
                    data_phase_0['IID'].append(sample)
                    data_phase_0[sample]['c_ancestry'] = lai.split(':')[1]
                    data_phase_0[sample]['start'] = pos

                # switch in ancestry --> store haplotype info --> chrom, start, end, ancestry, IID
                if data_phase_1[sample]['c_ancestry'] != lai.split(':')[2]:
                    data_phase_1['chrom'].append(chrom)
                    data_phase_1['start'].append(data_phase_1[sample]['start'] - 1)
                    data_phase_1['end'].append(data_phase_1[sample]['prev_pos'])
                    data_phase_1['la'].append(id_ancestry_mapping[data_phase_1[sample]['c_ancestry']])
                    data_phase_1['IID'].append(sample)
                    data_phase_1[sample]['c_ancestry'] = lai.split(':')[2]
                    data_phase_1[sample]['start'] = pos
                # update pointer to previous data point
                data_phase_0[sample]['prev_pos'] = pos
                data_phase_1[sample]['prev_pos'] = pos
    # save info about last haplotypes
    for sample, lai in zip(samples, line[9:]):
        if pos > data_phase_0[sample]['start'] + 1:
            data_phase_0['chrom'].append(chrom)
            data_phase_0['start'].append(data_phase_0[sample]['start'] - 1)
            data_phase_0['end'].append(data_phase_0[sample]['prev_pos'])
            data_phase_0['la'].append(id_ancestry_mapping[data_phase_0[sample]['c_ancestry']])
            data_phase_0['IID'].append(sample)
        if pos > data_phase_1[sample]['start'] + 1:
            data_phase_1['chrom'].append(chrom)
            data_phase_1['start'].append(data_phase_1[sample]['start'] - 1)
            data_phase_1['end'].append(data_phase_1[sample]['prev_pos'])
            data_phase_1['la'].append(id_ancestry_mapping[data_phase_1[sample]['c_ancestry']])
            data_phase_1['IID'].append(sample)

    # create df and save
    for sample in samples:
        del data_phase_0[sample]
        del data_phase_1[sample]
    df_phase_0 = pd.DataFrame.from_dict(data_phase_0)
    df_phase_1 = pd.DataFrame.from_dict(data_phase_1)
    df_phase_0.to_csv(args.output_phase0, sep='\t', header=True, index=False)
    df_phase_1.to_csv(args.output_phase1, sep='\t', header=True, index=False)


if __name__ == '__main__':
    main(sys.argv[1:])
