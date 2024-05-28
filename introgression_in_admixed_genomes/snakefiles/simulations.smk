def get_number_of_admixed_african_american_individuals(wildcards):
    n_indvs = 0
    with open(f'{data_path}AA_sample_ids.txt', 'r') as aa_ids:
        for line in aa_ids:
            n_indvs += 1
    aa_ids.close()
    return n_indvs

rule run_simulations:
    input:
        aa_ids = data_path + "AA_sample_ids.txt"
    output:
        vcfs=[simulations_path + f'chr{chrom}_replicate_' + '{n}.vcf.gz' for chrom in np.arange(1, 11)],
        trees=[simulations_path + '' + f'chr{chrom}_replicate_' + '{n}.trees' for chrom in np.arange(1, 11)],
        refs =simulations_path + 'reference_panel_replicate{n}.txt',
        afr_individuals = simulations_path + 'AFR_sample_ids_replicate{n}.txt',
        eur_individuals = simulations_path + 'EUR_sample_ids_replicate{n}.txt',
        eas_individuals = simulations_path + 'EAS_sample_ids_replicate{n}.txt',
        nea_individuals = simulations_path + 'neanderthal_sample_id_replicate{n}.txt',
        den_individuals = simulations_path + 'denisovan_sample_id_replicate{n}.txt',
        admixed_individuals = simulations_path + 'AA_sample_ids_replicate{n}.txt',
        genome_files =[simulations_path + '' + f"genomefile_chr{chrom}_replicate" + "{n}.bed" for chrom in np.arange(1, 11)],
        masks=[simulations_path + '' + f"regions_to_exclude_chr{chrom}_replicate" + "{n}.bed" for chrom in np.arange(1,11)]
    conda:
        '../envs/stdpopsim.yaml'
    params:
        vcf=lambda wildcards: simulations_path + 'chr{chrom}_replicate_' + f'{wildcards.n}.vcf.gz',
        trees=lambda wildcards: simulations_path + 'chr{chrom}_replicate_' + f'{wildcards.n}.trees',
        ref=simulations_path + 'reference_panel_replicate{n}.txt',
        afr_individuals=simulations_path + 'AFR_sample_ids_replicate{n}.txt',
        eur_individuals=simulations_path + 'EUR_sample_ids_replicate{n}.txt',
        eas_individuals=simulations_path + 'EAS_sample_ids_replicate{n}.txt',
        nea_individual=simulations_path + 'neanderthal_sample_id_replicate{n}.txt',
        den_individual=simulations_path + 'denisovan_sample_id_replicate{n}.txt',
        admixed_individuals=simulations_path + 'AA_sample_ids_replicate{n}.txt',
        genome_file=lambda wildcards: simulations_path + "genomefile_chr{chrom}_replicate" + f"{wildcards.n}.bed",
        mask_file=lambda wildcards: simulations_path + "regions_to_exclude_chr{chrom}_replicate" + f"{wildcards.n}.bed",
        n_amr=get_number_of_admixed_african_american_individuals,
    shell:
        'scripts/simulate_human_demography_with_archaic_introgression_and_american_admixture.py '
        '--trees {params.trees} --vcf {params.vcf} --ref {params.ref} '
        '-g {params.genome_file} --afr {params.afr_individuals} --eur {params.eur_individuals} '
        '--aa {params.admixed_individuals} --eas {params.eas_individuals} --nea {params.nea_individual} '
        '--den {params.den_individual} --n_amr {params.n_amr} --masks {params.mask_file}'

use rule window_genome as window_genomefile_sim with:
    input:
        simulations_path +  "genomefile_chr{chr}_replicate{n}.bed"
    output:
        simulations_path +  "chr{chr}_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed"
    params:
        bedtools = bedtools_path

rule combine_archaic_ids:
    input:
        nea_individual = simulations_path +  "neanderthal_sample_id_replicate{n}.txt",
        den_individual = simulations_path +  "denisovan_sample_id_replicate{n}.txt"
    output:
        simulations_path +  "archaic_sample_ids_replicate{n}.txt"
    shell:
        "cat {input} > {output}"

rule vcf_to_plink_sim:
    input:
        vcf=simulations_path + 'chr{chr}_replicate_{n}.vcf.gz',
        archaic=simulations_path +  "archaic_sample_ids_replicate{n}.txt"
    output:
        temp(multiext(simulations_path + 'tmp_chr{chr}_replicate_{n}', '.bed', '.bim', '.fam'))
    conda:
        '../envs/plink.yaml'
    params:
        prefix = simulations_path + 'tmp_chr{chr}_replicate_{n}'
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    threads: 8
    shell:
        'plink2 --vcf {input.vcf} --make-bed --set-all-var-ids @:# --max-alleles 2 --remove {input.archaic} --snps-only '
        '--maf 0.01 --out {params.prefix} --threads {threads} --memory {resources.mem_mb}'

rule get_new_ids_sim:
    input:
        fam=simulations_path + 'tmp_chr{chr}_replicate_{n}.fam',
        ref=simulations_path +  "reference_panel_replicate{n}.txt"
    output:
        temp(simulations_path +  "updated_fids_chr{chr}_replicate{n}.txt")
    shell:
        "grep -w -f <(cut -f2 {input.fam}) {input.ref} | "
        "awk 'BEGIN{{OFS=\"\\t\"}}{{print \"0\", $1, $2, $1}}' > {output}"

rule update_fids_sim:
    input:
        multiext(simulations_path + 'tmp_chr{chr}_replicate_{n}','.bed','.bim','.fam'),
        ids=simulations_path +  "updated_fids_chr{chr}_replicate{n}.txt"
    output:
        temp(multiext(simulations_path + 'chr{chr}_replicate_{n}','.bed','.bim','.fam'))
    params:
        input_prefix=simulations_path + 'tmp_chr{chr}_replicate_{n}',
        output_prefix=simulations_path + 'chr{chr}_replicate_{n}'
    conda:
        "../envs/plink.yaml"
    threads: 8
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    wildcard_constraints:
        n="[0-9]+"
    shell:
        "plink2 --bfile {params.input_prefix} --make-bed --out {params.output_prefix} --update-ids {input.ids} "
        "--threads {threads} --memory {resources.mem_mb}"

rule ld_prune_sim:
    input:
        multiext(simulations_path + 'chr{chr}_replicate_{n}','.bed','.bim','.fam')
    output:
        temp(multiext(simulations_path + 'chr{chr}_replicate_{n}','.prune.in','.prune.out'))
    params:
        input_prefix = simulations_path + 'chr{chr}_replicate_{n}'
    conda:
        "../envs/plink.yaml"
    threads: 16
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --bfile {params.input_prefix} --indep-pairwise 200 kb 1 0.1 --out {params.input_prefix} "
        "--threads {threads} --memory {resources.mem_mb}"
#
rule filter_ld_pruned_variants_sim:
    input:
        plink=multiext(simulations_path + 'chr{chr}_replicate_{n}','.bed','.bim','.fam'),
        keep=simulations_path + 'chr{chr}_replicate_{n}.prune.in'
    output:
        multiext(simulations_path + 'chr{chr}_replicate_{n}_ld_pruned','.bed','.bim','.fam')
    params:
        input_prefix = simulations_path + 'chr{chr}_replicate_{n}',
        output_prefix= simulations_path + 'chr{chr}_replicate_{n}_ld_pruned',
    conda:
        "../envs/plink.yaml"
    threads: 8
    wildcard_constraints:
        n="[0-9]+"
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        'plink2 --bfile {params.input_prefix} --extract {input.keep} --make-bed --out {params.output_prefix} '
        '--threads {threads} --memory {resources.mem_mb}'

def get_input_merge_chromosomes_replicates_sim(wildcards):
    bed  = [simulations_path +  f'chr{chrom}_replicate_{wildcards.n}_ld_pruned.bed' for chrom in np.arange(1, 11)]
    fam  = [simulations_path +  f'chr{chrom}_replicate_{wildcards.n}_ld_pruned.fam' for chrom in np.arange(1, 11)]
    bim  = [simulations_path +  f'chr{chrom}_replicate_{wildcards.n}_ld_pruned.bim' for chrom in np.arange(1, 11)]
    return {'bed': bed, "bim": bim, "fam": fam}

rule merge_chromosomes_replicates_sim:
    input:
        unpack(get_input_merge_chromosomes_replicates_sim)
    output:
        temp(multiext(simulations_path + 'ALL_chromosomes_replicate_{n}', ".bed", ".bim", ".fam")),
        bedfiles=temp("bed_files_replicate{n}.txt")
    params:
        output_base = simulations_path + 'ALL_chromosomes_replicate_{n}'
    threads: 8
    resources:
        mem_mb = 64000
    conda:
        "../envs/plink.yaml"
    wildcard_constraints:
        n="[0-9]+"
    shell:
        "for bed in {input.bed}; do echo $bed | sed 's/.bed//' >> {output.bedfiles}; done; "
        "plink --merge-list {output.bedfiles} --make-bed --allow-no-sex --out {params.output_base} --threads {threads} "
        "--memory {resources.mem_mb}"

rule compute_pca_reference_sim:
    input:
        multiext(simulations_path + 'ALL_chromosomes_replicate_{n}','.bed','.bim','.fam'),
    output:
        temp(multiext(simulations_path + 'ALL_chromosomes_replicate_{n}_ref',
             ".eigenvec", ".eigenval", '.acount', ".eigenvec.allele"))
    params:
        input_base = simulations_path + 'ALL_chromosomes_replicate_{n}',
        output_base = simulations_path + 'ALL_chromosomes_replicate_{n}_ref',
        pcs = M,
        maf = maf
    conda:
        "../envs/plink.yaml"
    threads: 16
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --bfile {params.input_base} --freq counts --pca {params.pcs} allele-wts --out {params.output_base} "
        "--threads {threads} --memory {resources.mem_mb} --maf {params.maf} "
        "--keep <(cat {params.input_base}.fam | cut -f1,2 | grep -v -E '^0')"

rule project_samples_sim:
    input:
        samples=multiext(simulations_path + 'ALL_chromosomes_replicate_{n}','.bed','.bim','.fam'),
        acount=simulations_path +  "ALL_chromosomes_replicate_{n}_ref.acount",
        wts=simulations_path +  "ALL_chromosomes_replicate_{n}_ref.eigenvec.allele"
    output:
        temp(simulations_path +  "ALL_chromosomes_replicate_{n}.sscore")
    conda:
        "../envs/plink.yaml"
    params:
        prefix=simulations_path + 'ALL_chromosomes_replicate_{n}',
        cols=6 + M - 1
    threads: 8
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --bfile {params.prefix} --read-freq {input.acount} --score {input.wts} 2 5 "
        "header-read no-mean-imputation variance-standardize --score-col-nums 6-{params.cols} --out {params.prefix} "
        "--threads {threads} --memory {resources.mem_mb}"

rule format_eigenvector_file_for_rye_sim:
    input:
        simulations_path +  "ALL_chromosomes_replicate_{n}.sscore"
    output:
        simulations_path +  "ALL_chromosomes_replicate_{n}.eigenvec"
    shell:
        "cat {input} | sed 's/_AVG//g' | cut -f1,2,5- > {output}"

rule copy_eigenval_sim:
    input:
        simulations_path +  "ALL_chromosomes_replicate_{n}_ref.eigenval"
    output:
        simulations_path +  "ALL_chromosomes_replicate_{n}.eigenval"
    shell:
        "cp {input} {output}"


rule create_rye_reference_sim:
    output:
        simulations_path +  "pop2group.txt"
    shell:
        "echo -e 'Pop\tSubgroup\tGroup' > {output}; "
        "echo -e 'EUR\tEuropean\tEuropean' >> {output}; "
        "echo -e 'AFR\tAfrican\tAfrican' >> {output}; "
        "echo -e 'EAS\tEastAsian\tEastAsian' >> {output}; "

rule determine_ancestry_proportions_sim:
    input:
        eigenvec=simulations_path +  "ALL_chromosomes_replicate_{n}.eigenvec",
        eigenval=simulations_path +  "ALL_chromosomes_replicate_{n}.eigenval",
        pop2group=simulations_path +  "pop2group.txt"
    output:
        simulations_path +  "ALL_chromosomes_replicate_{n}_ancestry_proportions" + f"-{M}.3.Q",
        simulations_path +  "ALL_chromosomes_replicate_{n}_ancestry_proportions" + f"-{M}.fam"
    params:
        rye = rye_path,
        out = simulations_path +  "ALL_chromosomes_replicate_{n}_ancestry_proportions",
        pcs = M
    threads: 16
    conda:
        "../envs/rye.yaml"
    shell:
        "{params.rye} --eigenval={input.eigenval} --eigenvec={input.eigenvec} --pop2group={input.pop2group} "
        "--output={params.out} --pcs={params.pcs} --threads={threads}"

rule generate_gt_and_run_ibdmix_simulations:
    input:
        vcf=simulations_path + 'chr{chr}_replicate_{n}.vcf.gz',
        archaic_indv=simulations_path +  "{archaic}_sample_id_replicate{n}.txt",
        archaic_all=simulations_path+ "archaic_sample_ids_replicate{n}.txt",
        afr=simulations_path +  "AFR_sample_ids_replicate{n}.txt",
        aa=simulations_path +  "AA_sample_ids_replicate{n}.txt",
        eas=simulations_path +  "EAS_sample_ids_replicate{n}.txt",
        eur=simulations_path +  "EUR_sample_ids_replicate{n}.txt",
        mask=simulations_path +  "regions_to_exclude_chr{chr}_replicate{n}.bed"
    output:
        eur=temp(simulations_path + "{archaic}_ibdmix_chr{chr}_EUR_results_replicate_{n}.tab"),
        afr=temp(simulations_path + "{archaic}_ibdmix_chr{chr}_AFR_results_replicate_{n}.tab"),
        eas=temp(simulations_path + "{archaic}_ibdmix_chr{chr}_EAS_results_replicate_{n}.tab"),
        aa=temp(simulations_path + "{archaic}_ibdmix_chr{chr}_AA_results_replicate_{n}.tab")
    params:
        ibdmix_dir = ibdmix_directory,
        bcftools = bcftools_path,
        archaic_error_rate= 0,
        max_modern_error_rate= 0,
        modern_error_proportion=0,
        min_minor_allele_count= min_minor_allele_count,
        lod_threshold= lod_threshold,
        gt=temp(simulations_path + '' +  "{archaic}_chr{chr}_ibdmix_genotypes_replicate_{n}"),
    resources:
        disk_mb = 220 * 1024
    wildcard_constraints:
        archaic = "neanderthal|denisovan",
        n="[0-9]+"
    shell:
        "{params.ibdmix_dir}generate_gt --modern <({params.bcftools} view -S ^{input.archaic_all} {input.vcf}) "
        "--archaic <({params.bcftools} view -S {input.archaic_indv} {input.vcf}) --output {params.gt}; "
        "{params.ibdmix_dir}ibdmix -g {params.gt} -o {output.eur} -s {input.eur} -r {input.mask} "
        "-m {params.min_minor_allele_count} -a {params.archaic_error_rate} -e {params.max_modern_error_rate} "
        "-c {params.modern_error_proportion} -d {params.lod_threshold}; "
        "{params.ibdmix_dir}ibdmix -g {params.gt} -o {output.afr} -s {input.afr} -r {input.mask} "
        "-m {params.min_minor_allele_count} -a {params.archaic_error_rate} -e {params.max_modern_error_rate} "
        "-c {params.modern_error_proportion} -d {params.lod_threshold}; "
        "{params.ibdmix_dir}ibdmix -g {params.gt} -o {output.eas} -s {input.eas} -r {input.mask} "
        "-m {params.min_minor_allele_count} -a {params.archaic_error_rate} -e {params.max_modern_error_rate} "
        "-c {params.modern_error_proportion} -d {params.lod_threshold}; "
        "{params.ibdmix_dir}ibdmix -g {params.gt} -o {output.aa} -s {input.aa} -r {input.mask} "
        "-m {params.min_minor_allele_count} -a {params.archaic_error_rate} -e {params.max_modern_error_rate} "
        "-c {params.modern_error_proportion} -d {params.lod_threshold}; rm {params.gt}"

# format IBDmix results
rule format_ibd_results_simulations:
    input:
        simulations_path +  "{archaic}_ibdmix_chr{chr}_{population}_results_replicate_{n}.tab"
    output:
        temp(simulations_path +  "{archaic}_ibdmix_chr{chr}_{population}_results_replicate_{n}_formatted.tab")
    wildcard_constraints:
        chr="[0-9]+",
    params:
        min_segment_length = minimum_length,
        lod_threshold = lod_threshold
    shell:
        "tail -n+2 {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($4 - $3 >= {params.min_segment_length} && "
        "$5 >= {params.lod_threshold}) print $2,$3,$4,$5,$1,\"{wildcards.population}\","
        "\"{wildcards.population}\"}}' > {output}"


def get_input_merge_ibd_results_simulations(wildcards):
    return [simulations_path +  f"{wildcards.archaic}_ibdmix_chr{wildcards.chr}_{population}_results_replicate_{wildcards.n}_formatted.tab"
            for population in ['AFR', 'EAS', 'AA', 'EUR']]

# merge IBDmix results for all samples
rule merge_ibd_results_simulations:
    input:
        get_input_merge_ibd_results_simulations,
    output:
        simulations_path + '{archaic}_introgressed_segments_chr{chr}_replicate_{n}.bed'
    params:
        min_segment_length = minimum_length,
        lod_threshold = lod_threshold,
        temp_dir=tmp_dir
    shell:
        "cat {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($3 - $2 >= {params.min_segment_length} && "
        "$4 >= {params.lod_threshold}) print $0}}' | sort -k1,1n -k2,2n -T {params.temp_dir} > {output}"

# mask Denisovan segments found in Africa in Neanderthal call set
use rule mask_denisovan_segments_in_neanderthal_set as mask_denisovan_segments_in_neanderthal_set_sim with:
    input:
        neanderthal = simulations_path + 'neanderthal_introgressed_segments_chr{chr}_replicate_{n}.bed',
        denisovan = simulations_path + 'denisovan_introgressed_segments_chr{chr}_replicate_{n}.bed'
    output:
        simulations_path + 'neanderthal_introgressed_segments_masked_denisovan_chr{chr}_replicate_{n}.bed'
    params:
        bedtools = bedtools_path,
        min_segment_length = minimum_length,
        lod_threshold = lod_threshold,
        afr_pop = "AFR"

def get_all_introgressed_segments(wildcards):
    return [simulations_path + f'neanderthal_introgressed_segments_masked_denisovan_chr{chrom}_replicate_{wildcards.n}.bed'
            for chrom in np.arange(1, 11)]

rule aggregate_nea_all:
    input:
        get_all_introgressed_segments
    output:
        simulations_path + 'neanderthal_introgressed_segments_masked_denisovan_replicate_{n}.bed'
    wildcard_constraints:
        n="[0-9]+",
    shell:
        'cat {input} | sort -k1,1n -k2,2n > {output}'

use rule mask_segments_in_AFR as mask_segments_in_AFR_sim with:
    input:
        simulations_path + 'neanderthal_introgressed_segments_masked_denisovan_replicate_{n}.bed'
    output:
        simulations_path + 'neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed'
    params:
        bedtools = bedtools_path

def get_genome_files(wildcards):
    return [simulations_path + f"genomefile_chr{chrom}_replicate{wildcards.n}.bed"
            for chrom in np.arange(1, 11)]

rule aggregate_genome_files:
    input:
        get_genome_files
    output:
        simulations_path + 'genomefile_replicate_{n}.bed'
    shell:
        "cat {input} | sort -k1,1n > {output}"

use rule window_genome as window_genomefile_sim_all with:
    input:
        simulations_path + 'genomefile_replicate_{n}.bed'
    output:
        simulations_path +  "genomefile_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed"
    params:
        bedtools = bedtools_path

rule calculate_introgressed_coverage_per_window_and_superpopulation_sim:
    input:
        sim=simulations_path +
                   "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}.bed",
        genomefile=simulations_path +  "genomefile_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed"
    output:
        temp(simulations_path +
             "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_coverage_per_window{windowsize}_s_{stepsize}.bed")
    params:
        bedtools=bedtools_path,
    wildcard_constraints:
        n="[0-9]+",
        windowsize="[0-9]+",
        stepsize="[0-9]+"
    shell:
        "{params.bedtools} intersect -a {input.genomefile} -b {input.sim} -wo | "
        "sort -k1,1n -k2,2n -k3,3n -k10 | {params.bedtools} groupby -g 1,2,3,10 -c 11 -o sum > {output}"

use rule calculate_introgressed_coverage_per_window_and_superpopulation_sim as calculate_introgressed_coverage_per_window_and_superpopulation_sim_masked with:
    input:
        sim=simulations_path +\
                   "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed",
        genomefile=simulations_path +  "genomefile_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed"
    output:
        temp(simulations_path +\
             "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed")
    params:
        bedtools=bedtools_path,
    wildcard_constraints:
        n="[0-9]+",
        windowsize="[0-9]+",
        stepsize="[0-9]+"

rule calculate_introgressed_coverage_per_window_and_individual_sim:
    input:
        sim=simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}.bed",
        genomefile=simulations_path +  "genomefile_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed"
    output:
        simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed"
    params:
        bedtools=bedtools_path,
    wildcard_constraints:
        windowsize="[0-9]+",
        stepsize="[0-9]+",
        n="[0-9]+"
    shell:
        "{params.bedtools} intersect -a {input.genomefile} -b {input.sim} -wo | "
        "sort -k1,1n -k2,2n -k3,3n -k8 | {params.bedtools} groupby -g 1,2,3,8 -c 11 -o sum > {output}"

use rule calculate_introgressed_coverage_per_window_and_individual_sim as calculate_introgressed_coverage_per_window_and_individual_sim_masked with:
    input:
        sim=simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed",
        genomefile=simulations_path +  "genomefile_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed"
    output:
        simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed"
    params:
        bedtools=bedtools_path,
    wildcard_constraints:
        windowsize="[0-9]+",
        stepsize="[0-9]+",
        n="[0-9]+"

# depending on demographic model not all super population show evidence of introgression
rule split_introgressed_coverage_per_window_and_superpopulation_sim:
    input:
        sim=simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_coverage_per_window{windowsize}_s_{stepsize}.bed"
    output:
        simulations_path + "{superpopulation}_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_coverage_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        superpopulation="AA|AFR|EUR|EAS",
        n="[0-9]+",
        windowsize="[0-9]+",
        stepsize="[0-9]+"
    shell:
        "if [[ $(grep -w {wildcards.superpopulation} {input.sim} | wc -l) == 0 ]]; "
        "then touch {output}; "
        "else grep -w {wildcards.superpopulation} {input.sim} > {output}; "
        "fi"

use rule split_introgressed_coverage_per_window_and_superpopulation_sim as split_introgressed_coverage_per_window_and_superpopulation_sim_masked with:
    input:
        sim=simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed"
    output:
        simulations_path + "{superpopulation}_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        superpopulation="AA|AFR|EUR|EAS",
        n="[0-9]+",
        windowsize="[0-9]+",
        stepsize="[0-9]+"

rule get_genetic_map_plink:
    input:
        data_path + genetic_map_path + "plink.formatted.chr16.GRCh38.map"
    output:
        simulations_path +  "plink.formatted.chr{chr}.GRCh38.map"
    params:
        chrom="{chr}"
    shell:
        "cat {input} | sed 's/chr16/chr{params.chrom}/' > {output}"

rule get_genetic_map:
    input:
        data_path + genetic_map_path + "genetic_map_Hg38_chr16.txt"
    output:
        simulations_path +  "genetic_map_Hg38_chr{chr}.txt"
    params:
        chrom="{chr}"
    shell:
        "cat {input} | sed 's/chr16/chr{params.chrom}/' > {output}"

rule train_flare_model_simulations:
    input:
        ref_vcf= simulations_path + 'chr1_replicate_{n}.vcf.gz',
        ref_panel= simulations_path +  "reference_panel_replicate{n}.txt",
        map= simulations_path +  "plink.formatted.chr1.GRCh38.map",
        target_vcf=simulations_path + 'chr1_replicate_{n}.vcf.gz',
        aa_individuals= simulations_path +  "AA_sample_ids_replicate{n}.txt",
    output:
        multiext(simulations_path + '' +  "african_american_and_ref_individuals_chr1_replicate{n}",
            ".model", ".anc.vcf.gz")
    params:
        base = simulations_path + '' +  "african_american_and_ref_individuals_chr1_replicate{n}",
        flare = flare_path,
        seed = 42,
        mem_gb= lambda wildcards,threads: threads * mem_gb_cpu
    threads: 16
    resources:
        mem_mb= lambda wildcards,threads: threads * mem_gb_cpu * 1000
    wildcard_constraints:
        n="[0-9]+",
    shell:
        "java -Xmx{params.mem_gb}g -jar {params.flare} ref={input.ref_vcf} ref-panel={input.ref_panel} "
        "gt={input.target_vcf} gt-samples={input.aa_individuals} "
        "map=<(sed 's/chr//' {input.map}) out={params.base} nthreads={threads} seed={params.seed}"

# paint chromosomes using previously trained flare model
rule run_flare_sim:
    input:
        ref_vcf = simulations_path + 'chr{chr}_replicate_{n}.vcf.gz',
        ref_panel = simulations_path +  "reference_panel_replicate{n}.txt",
        map = simulations_path +  "plink.formatted.chr{chr}.GRCh38.map",
        model = simulations_path + '' +  "african_american_and_ref_individuals_chr1_replicate{n}.model",
        target_vcf=simulations_path + 'chr{chr}_replicate_{n}.vcf.gz',
        aa_individuals= simulations_path +  "AA_sample_ids_replicate{n}.txt",
    output:
        simulations_path + '' +  "african_american_and_ref_individuals_chr{chr}_replicate{n}.anc.vcf.gz"
    wildcard_constraints:
        chr="[2-9][0-9]?|1[0-9]",
        n="[0-9]+",
    params:
        base = simulations_path + '' +  "african_american_and_ref_individuals_chr{chr}_replicate{n}",
        flare = flare_path,
        seed = 42,
        mem_gb = lambda wildcards,threads: threads * mem_gb_cpu
    threads: 16
    resources:
        mem_mb= lambda wildcards,threads: threads * mem_gb_cpu * 1000,
        load=10
    shell:
        "java -Xmx{params.mem_gb}g -jar {params.flare} ref={input.ref_vcf} ref-panel={input.ref_panel} "
        "gt={input.target_vcf} gt-samples={input.aa_individuals} "
        "map=<(sed 's/chr//' {input.map}) out={params.base} nthreads={threads} "
        "seed={params.seed} em=false model={input.model}"

use rule extract_local_ancestry as extract_local_ancestry_sim with:
    input:
        lai=simulations_path +  "african_american_and_ref_individuals_chr{chr}_replicate{n}.anc.vcf.gz",
    output:
        phase0=simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_phase0.bed',
        phase1=simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_phase1.bed'
    wildcard_constraints:
        n="[0-9]+",


rule compute_local_ancestry_coverage_per_window_and_phase_sim:
    input:
        genomefile=simulations_path +  "chr{chr}_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed",
        lai=simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_phase{phase}.bed'
    output:
        simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_per_window{windowsize}_s_{stepsize}_phase{phase}.bed'
    wildcard_constraints:
        n="[0-9]+",
        phase='0|1',
        windowsize="[0-9]+",
        stepsize="[0-9]+",
    params:
        bedtools = bedtools_path,
    shell:
        "{params.bedtools} intersect -a {input.genomefile} -b {input.lai} -wo | sort -k1,1n -k2,2n -k3,3n -k7 | "
        "{params.bedtools} groupby -g 1,2,3,7 -c 8 -o count > {output}"

rule compute_local_ancestry_frequencies_per_window_sim:
    input:
        lai_0=simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_per_window{windowsize}_s_{stepsize}_phase0.bed',
        lai_1=simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_per_window{windowsize}_s_{stepsize}_phase1.bed'
    output:
        simulations_path + 'chr{chr}_local_ancestry_information_replicate_{n}_per_window{windowsize}_s_{stepsize}.bed'
    wildcard_constraints:
        n="[0-9]+",
        windowsize="[0-9]+",
        stepsize="[0-9]+",
    params:
        bedtools = bedtools_path
    shell:
        "cat {input.lai_0} {input.lai_1} | sort -k1,1n -k2,2n -k3,3n -k4 | "
        "{params.bedtools} groupby -g 1-4 -c 5 -o sum  > {output}"

def get_lai_files(wildcards):
    return [simulations_path +
            f'chr{chrom}_local_ancestry_information_replicate_{wildcards.n}_per_window{wildcards.windowsize}_s_{wildcards.stepsize}.bed'
            for chrom in np.arange(1, 11)]

rule aggregate_lai_across_chromosomes:
    input:
        get_lai_files
    output:
        simulations_path + 'local_ancestry_information_replicate_{n}_per_window{windowsize}_s_{stepsize}.bed'
    shell:
        'cat {input} | sort -k1,1n -k2,2n > {output}'

def get_number_of_admixed_african_american_individuals_filtered(wildcards):
    n_indvs = 0
    with open(f'{simulations_path}AA_sample_ids_replicate{wildcards.n}.txt', 'r') as aa_ids:
        for line in aa_ids:
            n_indvs += 1
    aa_ids.close()
    return n_indvs

def get_masks_sim(wildcards):
    return [simulations_path + f"regions_to_exclude_chr{chrom}_replicate{wildcards.n}.bed"
            for chrom in np.arange(1, 11)]

def get_mask_pattern_sim(wildcards):
    return simulations_path + f"regions_to_exclude_chr" + "{chrom}" + f"_replicate{wildcards.n}.bed"

rule find_putatively_selected_neanderthal_segments_african_americans_sim:
    input:
        masks = get_masks_sim,
        lai=simulations_path + 'local_ancestry_information_replicate_{n}_per_window{windowsize}_s_{stepsize}.bed',
        nea_amr=simulations_path +  "AA_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed",
        nea_eur=simulations_path +  "EUR_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed",
        nea_eas=simulations_path +  "EAS_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed",
        introgressed_all=simulations_path + 'neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed',
        introgressed_all_coverage = simulations_path +  "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked" +
                                    "_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed",
        genomefile=simulations_path + 'genomefile_replicate_{n}.bed',
        windowed_genomefile=simulations_path + 'genomefile_windowed_w_{windowsize}_s_{stepsize}_replicate{n}.bed',
        aa_ids=simulations_path + "AA_sample_ids_replicate{n}.txt",
        genetic_maps=expand(simulations_path +  "genetic_map_Hg38_chr{chrom}.txt", chrom=np.arange(1, 11)),
        chrom_sizes=data_path + reference_path + 'hg38.chrom.sizes.bed'
    output:
        selected=simulations_path +  "AMR_putatively_selected_neanderthal_segments" +
                 "_replicate_{n}_window{windowsize}_s_{stepsize}.bed",
        all_pval=simulations_path +  "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked" +
                 "_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_pvalues.bed",
        all_exp=simulations_path +  "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked" +
                 "_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_expectations.bed"
    params:
        min_length = minimum_length_selected,
        max_eas_anc = 1 - min_afr_eur_component_combined,
        window_size = windowsize,
        step_size = stepsize,
        alpha = alpha,
        min_expectation = min_expectation,
        reps = 0,
        mask_pattern = get_mask_pattern_sim,
        tmp_dir= tmp_dir
    threads: 32
    wildcard_constraints:
        n="[0-9]+",
        windowsize=str(windowsize),
        stepsize=str(stepsize),
    shell:
        "scripts/identify_putatively_selected_introgressed_segments_in_african_americans.py -l {input.lai} "
        "--nea_amr {input.nea_amr} --nea_eur {input.nea_eur} --nea_eas {input.nea_eas} "
        "--ibdmix_all {input.introgressed_all} --ibdmix_all_coverage {input.introgressed_all_coverage} "
        "--min_length {params.min_length} --max_eas_anc {params.max_eas_anc} --min_expectation {params.min_expectation} "
        "-g {input.genomefile} -os {output.selected} -oip {output.all_pval} -oie {output.all_exp} "
        "-w {params.window_size} -s {params.step_size} --alpha {params.alpha} --reps {params.reps} --threads {threads} "
        "--genetic_maps {input.genetic_maps} --chrom_sizes {input.chrom_sizes} --mask_pattern {params.mask_pattern} "
        "--tmp_dir {params.tmp_dir} --windowed_genome_file {input.windowed_genomefile}"


def get_lai_pattern_sim(wildcards):
    return simulations_path + 'chr{chrom}_local_ancestry_information_' + \
           f'replicate_{wildcards.n}_' + 'phase{phase}.bed'

rule find_neanderthal_introgression_deserts_AMR_sim:
    input:
        masks=get_masks_sim,
        coverage=simulations_path + "AA_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_" +
                 f"coverage_per_window{windowsize}_s_{stepsize}.bed",
        coverage_masked=simulations_path + "AA_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_" +\
                        f"coverage_per_window{windowsize}_s_{stepsize}.bed",
        genomefile=simulations_path + 'genomefile_replicate_{n}.bed',
        ibdmix_all=simulations_path +
                   "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}.bed",
        ibdmix_masked=simulations_path +
                      "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed",
        lai=[simulations_path + '' + f'chr{chrom}_local_ancestry_information_' + 'replicate_{n}_' +
             f'phase{phase}.bed' for chrom in np.arange(1, 11)  for phase in [0, 1]],
        cw_eur= simulations_path + "EUR_introgression_frequencies_and_rank_callable_windows_replicate_{n}_afr_masked.bed",
        cw_eas= simulations_path + "EAS_introgression_frequencies_and_rank_callable_windows_replicate_{n}_afr_masked.bed"
    output:
        deserts=simulations_path + "AMR_introgression_deserts_replicate_{n}.bed",
        ranks=simulations_path + "AMR_introgression_frequencies_and_rank_callable_windows_replicate_{n}.bed",
        ranks_masked=simulations_path + "AMR_introgression_frequencies_and_rank_callable_windows_replicate_{n}_afr_masked.bed",
        deserts_new=simulations_path + "AMR_novel_introgression_deserts_replicate_{n}.bed",
        deserts_new_pvalues=simulations_path + "AMR_novel_introgression_deserts_pvalues_replicate_{n}.bed"
    params:
        mask_pattern = get_mask_pattern_sim,
        max_masked = max_masked,
        stride = stride,
        window_sizes = ' '.join([str(w) for w in window_sizes]),
        super_pop = 'AA',
        lai_pattern = get_lai_pattern_sim,
        tmp_dir = tmp_dir
    threads: 64
    resources:
        mem_mb=128 * 1000,
        load=5
    shell:
        'scripts/find_introgression_deserts.py --introgression_coverage {input.coverage} '
        '--introgression_coverage_masked {input.coverage_masked} --stride {params.stride} '
        '--max_masked {params.max_masked} --window_sizes {params.window_sizes} '
        '--genomefile {input.genomefile} --mask_pattern {params.mask_pattern} -o {output.deserts} '
        '--ibdmix {input.ibdmix_all} --ibdmix_masked {input.ibdmix_masked} --threads {threads} '
        '--output_ranks {output.ranks} --output_ranks_masked {output.ranks_masked} '
        '--super_pop {params.super_pop} --lai {params.lai_pattern} --callable_windows_references '
        '{input.cw_eur} {input.cw_eas} --callable_windows_references_labels EUR EAS '
        '--unique_deserts {output.deserts_new} --unique_deserts_pvalues {output.deserts_new_pvalues} '
        '--tmp_dir {params.tmp_dir}'

rule find_neanderthal_introgression_deserts_sim:
    input:
        masks=get_masks_sim,
        coverage=simulations_path + "{superpopulation}_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_" +
                 f"coverage_per_window{windowsize}_s_{stepsize}.bed",
        coverage_masked=simulations_path + "{superpopulation}_neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked_" +
                 f"coverage_per_window{windowsize}_s_{stepsize}.bed",
        genomefile=simulations_path + 'genomefile_replicate_{n}.bed',
        ibdmix_all=simulations_path +
                   "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}.bed",
        ibdmix_masked=simulations_path +
                      "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed",
    output:
        deserts=simulations_path + "{superpopulation}_introgression_deserts_replicate_{n}.bed",
        ranks=simulations_path + "{superpopulation}_introgression_frequencies_and_rank_callable_windows_replicate_{n}.bed",
        ranks_masked=simulations_path + "{superpopulation}_introgression_frequencies_and_rank_callable_windows_replicate_{n}_afr_masked.bed"
    params:
        mask_pattern = get_mask_pattern_sim,
        max_masked = max_masked,
        stride = stride,
        window_sizes = ' '.join([str(w) for w in window_sizes]),
        super_pop = '{superpopulation}',
        tmp_dir= tmp_dir
    wildcard_constraints:
        superpopulation = 'EUR|EAS'
    threads: 32
    resources:
        mem_mb=128 * 1000,
        load=5
    shell:
        'scripts/find_introgression_deserts.py --introgression_coverage {input.coverage} '
        '--introgression_coverage_masked {input.coverage_masked} --stride {params.stride} '
        '--max_masked {params.max_masked} --window_sizes {params.window_sizes} --output_ranks {output.ranks} '
        '--output_ranks_masked {output.ranks_masked} --genomefile {input.genomefile} '
        '--mask_pattern {params.mask_pattern} -o {output.deserts} --threads {threads} --super_pop {params.super_pop} '
        '--ibdmix {input.ibdmix_all} --ibdmix_masked {input.ibdmix_masked} --tmp_dir {params.tmp_dir}'
