#!/usr/bin/env python3
import tskit
import multiprocessing as mp
import os
import subprocess


def extraction_helper(args):
    chrom, rep = args
    if os.path.isfile(f'simulations/local/chr{chrom}_replicate_{rep}.vcf.gz'):
        return "done"
    else:
        ts = tskit.load(f'simulations/local/chr{chrom}_replicate_{rep}.trees')
        with open(f'simulations/local/chr{chrom}_replicate_{rep}.vcf', "w") as vcf_file:
            ts.write_vcf(vcf_file, contig_id=chrom)
        vcf_file.close()
        subprocess.run(["bgzip",  f'simulations/local/chr{chrom}_replicate_{rep}.vcf'])
        with open(f'simulations/local/genomefile_chr{chrom}_replicate{rep}.bed', 'w') as gf:
            gf.write(f"chr{chrom}\t{int(ts.sequence_length)}\n")
        gf.close()


def main():
    ready_to_map  = [(chrom, rep) for chrom in range(1, 11) for rep in range(10)]
    pool = mp.Pool(processes=10)
    pool.map(extraction_helper, ready_to_map)
    pool.close()
    pool.join()


if __name__ == '__main__':
    main()
