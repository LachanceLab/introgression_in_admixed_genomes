#!/usr/bin/env python3
import gzip


def write_line_and_prepare_for_next_segment(file, data_dict, sample, chrom, pos, ancestry, gt, alleles, phase):
    file.write(f'chr{chrom}\t{data_dict[sample][f"start_{phase}"]}\t'
                       f'{data_dict[sample][f"end_{phase}"]}\t{data_dict[sample][f"an{phase}"]}\t'
                       f'{sample}\t{",".join([str(site) for site in data_dict[sample][f"sites_{phase}"]])}\t'
                       f'{",".join(data_dict[sample][f"alleles_{phase}"])}\n')
    data_dict[sample][f'start_{phase}'] = pos - 1
    data_dict[sample][f'end_{phase}'] = pos
    data_dict[sample][f'an{phase}'] = ancestry
    data_dict[sample][f'alleles_{phase}'] = [alleles[gt]]
    data_dict[sample][f'sites_{phase}'] = [pos - 1]
    return data_dict


def update_data_dict(data_dict, pos, alleles, gt, phase):
    data_dict[sample][f'end_{phase}'] = pos
    data_dict[sample][f'alleles_{phase}'].append(alleles[gt])
    data_dict[sample][f'sites_{phase}'].append(pos - 1)
    return data_dict


ancestry_segments = dict()
phase_0_file = open('ancestry_segments_phase_0.bed', 'w')
phase_1_file = open('ancestry_segments_phase_1.bed', 'w')

with gzip.open('../ALL.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.anc.vcf.gz', 'rb') as vcf:
    for line in vcf:
        if line.startswith(b'##'):
            continue
        elif line.startswith(b'#CHROM'):
            samples = line.decode('utf-8').strip().split('\t')[9:]
            continue
        line = line.decode('utf-8').strip().split('\t')
        chrom = line[0]
        pos = int(line[1])
        ref = line[3]
        alt = line[4]
        alleles = [ref]
        alleles.extend(alt.split(','))
        for sample, format in zip(samples, line[9:]):
            phase_0 = int(format.split(":")[0].split('|')[0])
            phase_1 = int(format.split(":")[0].split('|')[1])
            an0 = int(format.split(":")[1])
            an1 = int(format.split(":")[2])

            if not sample in ancestry_segments:
                ancestry_segments[sample] = {'start_0': pos-1, 'end_0': pos, 'an0': an0,
                                             'alleles_0': [alleles[phase_0]], 'sites_0': [pos - 1],
                                             'start_1': pos - 1, 'end_1': pos, 'an1': an1,
                                             'alleles_1': [alleles[phase_1]], 'sites_1': [pos - 1]}
                continue

            if an0 == ancestry_segments[sample]['an0']:
                ancestry_segments = update_data_dict(ancestry_segments, pos, alleles, phase_0, 0)
            else:
                ancestry_segments = write_line_and_prepare_for_next_segment(phase_0_file, ancestry_segments, sample,
                                                                            chrom, pos, an0, phase_0, alleles, 0)

            if an1 == ancestry_segments[sample]['an1']:
                ancestry_segments = update_data_dict(ancestry_segments, pos, alleles, phase_1, 1)
            else:
                ancestry_segments = write_line_and_prepare_for_next_segment(phase_1_file, ancestry_segments, sample,
                                                                            chrom, pos, an1, phase_1, alleles, 1)

    for sample in samples:
        ancestry_segments = write_line_and_prepare_for_next_segment(phase_0_file, ancestry_segments, sample, chrom, pos,
                                                                    an0, phase_0, alleles, 0)
        ancestry_segments = write_line_and_prepare_for_next_segment(phase_1_file, ancestry_segments, sample, chrom, pos,
                                                                    an1, phase_1, alleles, 1)

vcf.close()
phase_0_file.close()
phase_1_file.close()
