import gzip
import re
from Bio import SeqIO

# create unique 35 bp reads:
checkpoint split_reference_into_35bp:
    input:
        seqbility=rules.download_and_install_seqbility.output,
        reference=rules.download_human_reference_genome.output
    output:
        directory(data_path + seqbility_tmp)
    shell:
        "mkdir -p {output}; {input.seqbility}/splitfa {input.reference} 35 | split -l 20000000 - {output}/split_"

# bwa index human reference genome
rule bwa_index_reference:
    input:
        rules.download_human_reference_genome.output
    output:
        multiext(data_path + reference_path + hg38_fasta_url.split('/')[-1].replace('.gz', ''),
            '.bwt', '.pac', '.ann', '.amb', '.sa')
    params:
        bwa = bwa_path
    shell:
        "{params.bwa} index -a bwtsw {input}"

# align all 35 kmers to human reference genome
rule bwa_align:
    input:
        split=data_path + seqbility_tmp + 'split_{part}',
        reference=rules.download_human_reference_genome.output,
        index=rules.bwa_index_reference.output
    output:
        temp(data_path + seqbility_tmp + 'mapped_{part}.sai')
    params:
        bwa = bwa_path
    resources:
        load=5
    shell:
        "{params.bwa} aln -R 1000000 -O 3 -E 3 {input.reference} {input.split} > {output}"

# convert alignments to sam files
rule bwa_samse:
    input:
        sai=data_path + seqbility_tmp + 'mapped_{part}.sai',
        split=data_path + seqbility_tmp + 'split_{part}',
        reference=rules.download_human_reference_genome.output
    output:
        temp(data_path + seqbility_tmp + 'mapped_{part}.sam')
    resources:
        load = 5
    params:
        bwa = bwa_path
    shell:
        "{params.bwa} samse {input.reference} {input.sai} {input.split} > {output}"

# get all sam files
def generate_raw_mappability_mask_input(wildcards):
    checkpoints.split_reference_into_35bp.get(**wildcards)
    parts = glob_wildcards(data_path + seqbility_tmp + 'split_{part}').part
    return expand(data_path + seqbility_tmp + 'mapped_{part}.sam', part=sorted(parts))

# generate raw mappability mask
rule generate_raw_mappability_mask:
    input:
        sams=generate_raw_mappability_mask_input,
        seqbility=rules.download_and_install_seqbility.output
    output:
        temp(data_path + mask_path + "human_genome_raw_mappability_mask_35.fa")
    params:
        data_path = data_path,
        seqbility_tmp=data_path + seqbility_tmp
    shell:
        "cat {input.sams} | {input.seqbility}/gen_raw_mask.pl > {output}; rm -r {params.seqbility_tmp}"

# at at least 50% of all possible 35mers overlapping a position do not find a match to any other position in the genome
# allowing for up to one mismatch
rule generate_mappability_mask:
    input:
        seqbility=rules.download_and_install_seqbility.output,
        raw_mask=rules.generate_raw_mappability_mask.output
    output:
        temp(data_path + mask_path + "human_genome_mappability_mask_35_50.fa")
    shell:
        "{input.seqbility}/gen_mask -l 35 -r 0.5 {input.raw_mask} > {output}"

# convert mappability mask from fasta format to bed format
rule mappability_mask_fa_to_bed:
    input:
        data_path + mask_path + "human_genome_mappability_mask_35_50.fa"
    output:
        data_path + mask_path + "human_genome_mappability_mask_35_50.bed"
    run:
        # parse fasta file
        mask_fa = SeqIO.to_dict(SeqIO.parse(input[0], 'fasta'))
        mask_bed = open(output[0], 'w+')
        for chrom in chromosomes:
            # get seq
            seq = mask_fa[f'chr{chrom}'].seq
            mappable = False
            for i, base in enumerate(seq):
                # start of mappable regions
                if not mappable and base == '3':
                    region_start = i
                    mappable = True
                # end of mappable region
                elif mappable and base != '3':
                    mask_bed.write(f'chr{chrom}\t{region_start}\t{i}\n')
                    mappable = False
            # if end of chromosome is mappable
            if mappable:
                mask_bed.write(f'chr{chrom}\t{region_start}\t{i + 1}\n')
        mask_bed.close()

# extract positions where humans have a derived allele
rule extract_snps_from_whole_genome_alignments:
    input:
        data_path + whole_genome_alignments_path + whole_genome_alignments_url.split('/')[-1]
    output:
        temp(data_path + whole_genome_alignments_path + "snps_human_{species}.bed")
        # data_path+ whole_genome_alignments_path + "chr{chr}.snps_human_{species}.bed"
    run:
        # open files
        with gzip.open(input[0], 'rt') as axtfile, open(output[0], 'wt') as bedfile:
            # iterate over all alignments
            while True:
                # parse, look https://genome.ucsc.edu/goldenPath/help/axt.html for format description
                summary = axtfile.readline()
                if summary.startswith('#'):
                    continue
                if summary == '':
                    return
                chrom = summary.split(' ')[1]
                start = int(summary.split(' ')[2]) - 1
                # unplaced chromosome
                if not re.match('^chr[0-9][0-9]?$', chrom):
                    axtfile.readline()
                    axtfile.readline()
                    axtfile.readline()
                else:
                    # read alignments
                    human = axtfile.readline().strip().upper()
                    other = axtfile.readline().strip().upper()
                    axtfile.readline()
                    # chrom = chrom.replace('chr', '')
                    gaps = 0
                    for i, (hum, oth) in enumerate(zip(human, other)):
                        # gaps in alignment --> skip
                        if hum == '-':
                            gaps += 1
                            continue
                        # humans have a different allele
                        if hum != oth and oth != '-':
                            bedfile.write(f"chr{chrom}\t{start + i - gaps}\t{hum}\t{oth}\n")
        axtfile.close()
        bedfile.close()

# get positions at which Africans have a C/G alternative allele, positions are 0 based --> $2 - 1
rule get_modern_human_cg_snps:
    input:
        vcf = data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        # vcf = data_path + merged_datasets_path + "1KG_HGDP_chr{chr}.vcf.gz",
        panel = data_path + "1000G_phase3/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        temp(data_path + mask_path + "chr{chr}.cg_mutations_modern.bed")
    params:
        bcftools = bcftools_path,
        afr_pop= "ESN|GWD|LWK|MSL|YRI"
    shell:
        "grep -E \"{params.afr_pop}\" {input.panel} | cut -f1 | {params.bcftools} view -S - -v snps -H {input.vcf} | "
        "cut -f1,2,4,5 | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{split($4, array, \",\")}}"
        "{{for (i in array) if (array[i] ~/C|G/ && length(array[i]) == 1 && length($3) == 1)  print $1,$2 - 1,$3,array[i]}}' "
        "> {output}"

# generate CpG mask, i.e., if C or G in one genome (human, other primates) is next to a C or G
# in another genome (human, other primates)
rule generate_cpg_mask:
    input:
        primates = expand(data_path + whole_genome_alignments_path + "snps_human_{species}.bed",species=species),
        # primates = expand(data_path + whole_genome_alignments_path + "chr{chr}.snps_human_{species}.bed", species=species, chr=chromosomes),
        modern = expand(data_path + mask_path + "chr{chr}.cg_mutations_modern.bed", chr=chromosomes),
        # modern = 'snps_human_modern.bed',
        reference = data_path + reference_path + hg38_fasta_url.split('/')[-1].replace('.gz', '')
    output:
        temp(data_path + mask_path + 'tmp_cpg.bed')
    shell:
        "scripts/generate_cpg_mask.py -i {input.primates} {input.modern} -r {input.reference} -o {output}"

# merge overlapping intervals in CpG mask
rule merge_intervals_in_cpg_mask:
    input:
        data_path + mask_path + 'tmp_cpg.bed'
    output:
        data_path + mask_path + 'cpg.bed'
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} merge -i {input} > {output}"

# remove all sites within 5bp of indels in reference files
rule create_indel_mask_1kgp:
    input:
        vcf=data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        pops= data_path + "1000G_phase3/individuals_in_reference_populations_1kgp.txt"
    output:
        temp(data_path + mask_path + "chr{chr}_indel_mask_1kgp.bed")
    params:
        bcftools = bcftools_path
    shell:
        # if deletion add the length of reference allele to window, if insertion just add window
        "{params.bcftools} view -v indels -H -S {input.pops} {input.vcf} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if (length($4) > 1) "
        "print $1, $2 - 6, $2 + 5 + length($4); else if (length($5) > 1) print $1, $2 - 6, $2 + 5}}'"
        " > {output}"

# use rule create_indel_mask_1kgp as create_indel_mask_hgdp with:
#     input:
#         data_path + "hgdp/hgdp_wgs.20190516.native_americans.chr{chr}.vcf.gz",
#     output:
#         temp(data_path + mask_path + "chr{chr}_indel_mask_hgdp.bed")
#     params:
#         bcftools = bcftools_path

rule merge_indel_mask:
    input:
        # data_path + mask_path + "chr{chr}_indel_mask_hgdp.bed",
        data_path+ mask_path + "chr{chr}_indel_mask_1kgp.bed"
    output:
        temp(data_path + mask_path + "chr{chr}_indel_mask.bed")
    params:
        bedtools = bedtools_path
    shell:
        "bedtools merge -i <(cat {input} | sort -k1,1 -k2,2n) > {output}"


# intersect accessibility masks for 1KGP and archaic genomes
rule intersect_accessible_masks:
    input:
        human = rules.lift_1kgp_mask_from_hg19_to_hg38.output,
        archaic = data_path + "archaic_genomes/{archaic_genome}/filter_bed/chr{chr}_mask.bed.gz",
        # archaic= data_path + "archaic_genomes/{archaic_genome}/filter_bed/{archaic_genome}.map35_50.MQ30.Cov.indels.TRF.bed",
        mappable = data_path + mask_path + "human_genome_mappability_mask_35_50.bed"
    output:
        temp(data_path + mask_path + "chr{chr}.accessible_and_mappable_mask.{archaic_genome}.bed")
        # temp(data_path + mask_path + "accessible_and_mappable_mask.{archaic_genome}.bed")
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} intersect -a {input.human} -b {input.archaic} | {params.bedtools} intersect -a - -b {input.mappable} "
        "> {output}"

# subtract accessible mask from entire genome
rule subtract_accessible_mask_from_entire_genome:
    input:
        chrom_sizes = rules.download_chromosome_sizes.output,
        incl_mask = data_path + mask_path + "chr{chr}.accessible_and_mappable_mask.{archaic_genome}.bed"
        # incl_mask= data_path + mask_path + "accessible_and_mappable_mask.{archaic_genome}.bed"
    output:
        temp(data_path + mask_path + "chr{chr}.non_accessible_or_mappable_mask.human.{archaic_genome}.bed")
    params:
        bedtools = bedtools_path
    shell:
        "cat {input.chrom_sizes} | grep -w \"^chr{wildcards.chr}\" | {params.bedtools} subtract -a - -b {input.incl_mask} "
        "> {output}"

# merge masks with regions to exclude (segmental duplications, 5bp within indels in modern, cpg islands, repeats)
rule merge_masks_to_exclude:
    input:
        dup = rules.download_segmental_duplication_map.output,
        indels = data_path + mask_path + "chr{chr}_indel_mask.bed",
        cpg = data_path + mask_path + "cpg.bed",
        repeats = rules.download_simple_repeat_map.output,
        not_callable = data_path + mask_path + "chr{chr}.non_accessible_or_mappable_mask.human.{archaic_genome}.bed",
        human_introgressed = data_path + mask_path + "human_to_{archaic_genome}_introgressed_regions.bed",
        gaps=rules.download_gaps.output
    output:
        data_path + mask_path + "chr{chr}.regions_to_exclude.{archaic_genome}.bed"
    params:
        bedtools = bedtools_path
    shell:
        "cat {input} | grep -w '^chr{wildcards.chr}' | sed 's/^chr//' | sort -k1,1n -k2,2n | {params.bedtools} merge -i - > {output}"