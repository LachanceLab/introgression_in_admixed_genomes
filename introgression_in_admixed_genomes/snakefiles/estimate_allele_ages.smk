def get_archaic_genomes(wildcards):
    return [data_path + f"archaic_genomes/{archaic_genome}/vcf/chr{wildcards.chrom}_mq25_mapab100.vcf.gz" if archaic_genome != 'Chagyrskaya'
            else arg_dir + f'{archaic_genome}_chr{wildcards.chrom}_mq25_mapab100.vcf.gz'
            for archaic_genome in archaic_genomes]


def get_command_archaic_genomes(wildcards):
    archaics = ''
    for archaic_genome in archaic_genomes:
        if archaic_genome == 'Altai':
            archaics += f"--vcf_altai {data_path}archaic_genomes/Altai/vcf/chr{wildcards.chrom}_mq25_mapab100.vcf.gz "
        if archaic_genome == 'Vindija33.19':
            archaics += f"--vcf_vindija {data_path}archaic_genomes/Vindija33.19/vcf/chr{wildcards.chrom}_mq25_mapab100.vcf.gz "
        if archaic_genome == 'Chagyrskaya':
            archaics += f"--vcf_chagyrskaya {arg_dir}Chagyrskaya_chr{wildcards.chrom}_mq25_mapab100.vcf.gz "
        if archaic_genome == 'Denisova':
            archaics += f"--vcf_denisova {data_path}archaic_genomes/Denisova/vcf/chr{wildcards.chrom}_mq25_mapab100.vcf.gz "
    return archaics


rule bgzip_chagyrskaya:
    input:
        data_path + "archaic_genomes/Chagyrskaya/vcf/chr{chrom}_mq25_mapab100.vcf.gz"
    output:
        vcf=temp(arg_dir + "Chagyrskaya_chr{chrom}_mq25_mapab100.vcf.gz"),
        tbi=temp(arg_dir + "Chagyrskaya_chr{chrom}_mq25_mapab100.vcf.gz.tbi")
    params:
        bcftools = bcftools_path
    threads: 8
    conda:
        "../envs/haplotype_dating.yaml"
    shell:
        "{params.bcftools} view -Oz -o {output.vcf} --threads {threads} {input}; tabix {output.vcf}"

checkpoint parse_vcfs_to_sample_data:
    input:
        get_archaic_genomes,
        vcf_modern=data_path + "1000G_phase3/REF.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_annotated.vcf.gz",
        tbi_modern=data_path + "1000G_phase3/REF.chr{chrom}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_annotated.vcf.gz.tbi",
        afr_sample_ids=data_path + 'AFR_sample_ids.txt',
        eur_sample_ids=data_path + 'EUR_sample_ids.txt',
        eas_sample_ids=data_path + 'EAS_sample_ids.txt'
    output:
        temp(directory(arg_dir + "sampledata_{chrom}"))
    params:
        generation_time = generation_time,
        output_dir = arg_dir + "sampledata_{chrom}/",
        archaics=get_command_archaic_genomes
    conda:
        "../envs/haplotype_dating.yaml"
    threads: 16
    shell:
        'scripts/estimate_allele_ages.py parse_vcfs --vcf_modern {input.vcf_modern} {params.archaics} '
        '--afr_ids {input.afr_sample_ids} --eur_ids {input.eur_sample_ids} --eas_ids {input.eas_sample_ids} '
        '-o {params.output_dir} --generation_time {params.generation_time} --chromosome {wildcards.chrom} '
        '--threads {threads}'


def get_input_merge_sample_data(wildcards):
    checkpoints.parse_vcfs_to_sample_data.get(**wildcards)
    parts = glob_wildcards(arg_dir + 'sampledata_{chrom}/partial_file_part{part}.samples').part
    return expand(arg_dir + 'sampledata_{chrom}/partial_file_part{part}.samples', part=sorted(list(set(parts))),
                  chrom=wildcards.chrom)

rule merge_sample_data:
    input:
        get_input_merge_sample_data
    output:
        temp(arg_dir + 'chr{chrom}.samples')
    params:
        output_dir = arg_dir
    conda:
        '../envs/haplotype_dating.yaml'
    wildcard_constraints:
        chrom="[0-9]+"
    shell:
        "scripts/estimate_allele_ages.py merge_sampledata --sampledata {input} "
        "--output_dir {params.output_dir} --chromosome {wildcards.chrom}"

rule infer_chromosome_arm:
    input:
        samples=arg_dir + 'chr{chrom}.samples',
        hapmap=data_path + genetic_map_path + "genetic_map_Hg38_chr{chrom}.txt",
        chromosome_arms=data_path + reference_path + 'hg38_chromosome_arms.bed'
    output:
        temp(arg_dir + 'chr{chrom}_{arm}_arm.trees'),
        temp(arg_dir+ 'chr{chrom}_{arm}_arm.samples')
    conda:
        '../envs/haplotype_dating.yaml'
    params:
        output_dir = arg_dir
    threads: 16
    shell:
        "scripts/estimate_allele_ages.py infer_chromosome_arm --sampledata {input.samples} "
        "--hapmap {input.hapmap} --chromosome_arms_df {input.chromosome_arms} --chromosome_arm {wildcards.arm} "
        "--output_dir {params.output_dir} --chromosome {wildcards.chrom} --threads {threads}"

rule date_alleles:
    input:
        arg_dir + 'chr{chrom}_{arm}_arm.trees'
    output:
        temp(arg_dir + 'chr{chrom}_{arm}_arm_dated.trees')
    conda:
        '../envs/haplotype_dating.yaml'
    params:
        mu=mutation_rate,
        Ne=effective_population_size,
    threads: 16
    shell:
        "scripts/estimate_allele_ages.py date_alleles --trees {input} --ne {params.Ne} -m {params.mu} "
        "--threads {threads}"

rule get_dated_sampledata:
    input:
        samples=arg_dir+ 'chr{chrom}_{arm}_arm.samples',
        trees=arg_dir + 'chr{chrom}_{arm}_arm_dated.trees'
    output:
        arg_dir + 'chr{chrom}_{arm}_arm_dated.samples'
    conda:
        '../envs/haplotype_dating.yaml'
    shell:
        "scripts/estimate_allele_ages.py get_dated_sampledata_from_ts --sampledata {input.samples} "
        "--trees {input.trees}"

rule generate_ancestors_ts:
    input:
        samples=arg_dir + 'chr{chrom}_{arm}_arm_dated.samples',
        hapmap=data_path + genetic_map_path + "genetic_map_Hg38_chr{chrom}.txt",
        chromosome_arms=data_path + reference_path + 'hg38_chromosome_arms.bed'
    output:
        temp(arg_dir + 'chr{chrom}_{arm}_arm_ancestors.trees')
    conda:
        '../envs/haplotype_dating.yaml'
    params:
        output_dir = arg_dir
    threads: 16
    shell:
        "scripts/estimate_allele_ages.py generate_ancestors_ts --sampledata {input.samples} "
        "--chromosome_arms_df {input.chromosome_arms} --hapmap {input.hapmap} --chromosome {wildcards.chrom} "
        "--chromosome_arm {wildcards.arm} --threads {threads} -o {params.output_dir}"

rule reinfer_ts:
    input:
        samples=arg_dir + 'chr{chrom}_{arm}_arm_dated.samples',
        trees=arg_dir + 'chr{chrom}_{arm}_arm_ancestors.trees',
        hapmap=data_path + genetic_map_path + "genetic_map_Hg38_chr{chrom}.txt",
        chromosome_arms=data_path + reference_path + 'hg38_chromosome_arms.bed'
    output:
        arg_dir + 'chr{chrom}_{arm}_arm_reinferred.trees'
    conda:
        '../envs/haplotype_dating.yaml'
    params:
        output_dir = arg_dir
    threads: 16
    shell:
        "scripts/estimate_allele_ages.py reinfer_ts --sampledata {input.samples} --trees {input.trees} "
        "--chromosome_arms_df {input.chromosome_arms} --hapmap {input.hapmap} --chromosome {wildcards.chrom} "
        "--chromosome_arm {wildcards.arm} --threads {threads} -o {params.output_dir}"

rule get_allele_ages:
    input:
        samples=arg_dir + 'chr{chrom}_{arm}_arm_dated.samples',
        trees=arg_dir + 'chr{chrom}_{arm}_arm_reinferred.trees'
    output:
        temp(arg_dir + 'chr{chrom}_{arm}_arm_tmrca_summarized.tab')
    conda:
        '../envs/haplotype_dating.yaml'
    params:
        generation_time=generation_time
    shell:
        "scripts/estimate_allele_ages.py get_allele_ages --sampledata {input.samples} --trees {input.trees} "
        "--chromosome {wildcards.chrom} -g {params.generation_time}"

def get_input_merge_allele_ages_per_chrom(wildcards):
    if wildcards.chrom != '13' and wildcards.chrom != '14' and wildcards.chrom != '15' and wildcards.chrom != '22':
        return [arg_dir + f'chr{wildcards.chrom}_p_arm_tmrca_summarized.tab',
                arg_dir + f'chr{wildcards.chrom}_q_arm_tmrca_summarized.tab']
    else:
        # no variants on p arm of chromosomes 13,14,15 and 22 in 1KGP data
        return [arg_dir + f'chr{wildcards.chrom}_q_arm_tmrca_summarized.tab']

rule merge_allele_ages_per_chrom:
    input:
        get_input_merge_allele_ages_per_chrom
        # p_arm=arg_dir + 'chr{chrom}_p_arm_tmrca_summarized.tab',
        # q_arm=arg_dir + 'chr{chrom}_q_arm_tmrca_summarized.tab'
    output:
        arg_dir + 'chr{chrom}_tmrca_summarized.tab'
    shell:
        "cat {input} | sort -k1,1 -k2,2n > {output}"
