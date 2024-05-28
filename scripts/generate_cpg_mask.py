#!/usr/bin/env python3
# Identify CpG sites by considering the human reference genome (hg38), human genomes, other primate genomes
# CpG sites are sites at which a site with a C or G in one genome is adjacent to a site with a C or G in another genome
# This approach was also used by Vernot et al. 2014 https://www.science.org/doi/10.1126/science.1245938
# The script is contains elements from a script by Chen et al. https://github.com/PrincetonUniversity/IBDmix/blob/main/j.cell.2020.01.012-workflow/generate_cpg_mask.py

import sys
import argparse
from Bio import SeqIO
import numpy as np
from collections import defaultdict


class Snps:
    def __init__(self, infile):
        self.data = defaultdict(dict)
        with open(infile, 'r') as bed:
            for line in bed:
                chrom, pos, ref, alt = line.strip().split('\t')
                pos = int(pos)
                self.data[chrom][pos] = (ref, alt)

    def is_cpg(self, chrom, position, base):
        try:
            ref, alt = self.data[chrom][position]
            assert (ref == base or ref in ('N', 'M') or base in ('N', 'M')), f'Expected {ref} at {chrom} {position} ' \
                                                                             f'to match {base}'
            return alt == 'C', alt == 'G'
        # site is not present
        except KeyError:
            return False, False


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_files', nargs='+',
                        help='Input files with SNPs where selected primate species and/or modern humans are '
                             'polymorphic. Expected format is:CHROM\tPOS\tREF\tALT. Expected file name is:'
                             '"[chr{chr}.cg|snps]_[mutations|human]_[{species}|modern].bed"')
    parser.add_argument('-r', '--reference', help='Human reference genome in FASTA format. [./data/reference/hg19.fa',
                        default='./data/reference/hg19.fa')
    parser.add_argument('-o', '--output', help='Output file to write identified CpG nucleotides to'
                                               '[./data/masks/cpg.bed]', default='./data/masks/cpg.bed')
    args = parser.parse_args()
    reference = SeqIO.to_dict(SeqIO.parse(args.reference, 'fasta'))
    chromosomes = np.arange(1, 23)
    snps = dict()
    outputfile = open(args.output, 'w')
    for chrom in chromosomes:
        for infile in args.input_files:
            if infile.split('/')[-1].startswith(f'chr{chrom}.'):
                 species = infile.split('/')[-1].split('_')[2].split('.bed')[0]
                 snps[species] = Snps(infile)
            elif chrom == 1 and not infile.split('/')[-1].startswith('chr'):
                 species = infile.split('/')[-1].split('_')[2].split('.bed')[0]
                 snps[species] = Snps(infile)
        sequence = reference[f'chr{chrom}'].seq.upper()
        last_base_C = False
        last_base_G = False
        for position, base in enumerate(sequence):
            base_is_C = base == 'C'
            base_is_G = base == 'G'
            for species, snp in snps.items():
                C, G = snp.is_cpg(str(chrom), position, base)
                base_is_C |= C
                base_is_G |= G
            # there is C in any genome and an adjacent G in any genome or
            # there is G in any genome and an adjacent C in any genome
            if (last_base_C and base_is_G) or (last_base_G and base_is_C):
                outputfile.write(f'chr{chrom}\t{position - 1}\t{position + 1}\n')
            last_base_C = base_is_C
            last_base_G = base_is_G
    outputfile.close()


if __name__ == '__main__':
    main(sys.argv[1:])

