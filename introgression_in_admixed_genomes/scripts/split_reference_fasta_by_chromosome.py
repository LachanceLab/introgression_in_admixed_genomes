#!/usr/bin/env python3
import sys
import numpy as np
import argparse
from Bio import SeqIO


def main(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', help='Input fasta')
    parser.add_argument('--prefix', help="Output prefix. Output file name will be {prefix}chr{chrom}{suffix}")
    parser.add_argument('--suffix', help="Output suffix. Output file name will be {prefix}chr{chrom}{suffix}")
    args = parser.parse_args()
    record_dict = SeqIO.to_dict(SeqIO.parse(args.input, "fasta"))
    for chrom in np.arange(1, 23):
        with open(f'{args.prefix}chr{chrom}{args.suffix}', "w") as output_handle:
            SeqIO.write(record_dict[f'chr{chrom}'], output_handle, 'fasta')
        output_handle.close()


if __name__ == "__main__":
    main(sys.argv[1:])