from Bio import SeqIO

# Generate genotype
rule generate_gt_ibdmix_1kgp:
    input:
        modern_vcf= data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        archaic_vcf = data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz"
    output:
        temp(results_path + ibdmix_genotypes_path + "chr{chr}_{archaic_genome}_ibdmix_genotypes_1kgp")
    params:
        ibdmix_dir = ibdmix_directory
    shell:
        "{params.ibdmix_dir}generate_gt --modern <(zcat {input.modern_vcf} | sed 's/^chr//') --archaic <(zcat {input.archaic_vcf} | sed 's/^chr//') "
        "--output {output}"


use rule generate_gt_ibdmix_1kgp as generate_gt_ibdmix_target with:
    input:
        modern_vcf= data_path + 'target_individuals_chr{chr}.vcf.gz',
        archaic_vcf = data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz"
    output:
        temp(results_path + ibdmix_genotypes_path + "chr{chr}_{archaic_genome}_ibdmix_genotypes_target")
    params:
        ibdmix_dir = ibdmix_directory


rule generate_gt_ibdmix_aou:
    input:
        modern_vcf = GS.remote(GS_PREFIX + '/' + all_of_us_vcf),
        tbi = GS.remote(GS_PREFIX + '/' + all_of_us_vcf + ".tbi"),
        archaic_vcf= data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz",
    output:
        geno=temp(results_path + ibdmix_genotypes_path + "chr{chr}_{archaic_genome}_ibdmix_genotypes_aou")
    params:
        ibdmix_dir=ibdmix_directory
    shell:
        "{params.ibdmix_dir}generate_gt --modern <(zcat {input.modern_vcf} | sed 's/^chr//') --archaic <(zcat {input.archaic_vcf} | sed 's/^chr//') "
        "--output {output.geno}; rm {input.modern_vcf}; rm {input.tbi}"




# use rule generate_gt_ibdmix_reference as generate_gt_ibdmix_target with:
#     input:
#         modern_vcf= data_path + 'target_individuals_chr{chr}.vcf.gz',
#         archaic_vcf = data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz"
#     output:
#         results_path + ibdmix_genotypes_path + "target_chr{chr}_{archaic_genome}_ibdmix_genotypes"
#     params:
#         ibdmix_dir = ibdmix_directory

# def get_gt_run_ibdmix(wildcards):
#     if wildcards.population == 'AA':
#         return results_path + ibdmix_genotypes_path + f"target_chr{wildcards.chr}_{wildcards.archaic_genome}_ibdmix_genotypes"
#     else:
#         return results_path + ibdmix_genotypes_path + f"reference_chr{wildcards.chr}_{wildcards.archaic_genome}_ibdmix_genotypes"

# detect introgressed segments using IBDmix
rule run_ibdmix_1kgp:
    input:
        gt = results_path + ibdmix_genotypes_path + "chr{chr}_{archaic_genome}_ibdmix_genotypes_1kgp",
        samples = data_path + "{population}_sample_ids.txt",
        mask = data_path + mask_path + "chr{chr}.regions_to_exclude.{archaic_genome}.bed",
        kin=related_samples
    output:
        temp(results_path + "ibdmix_{archaic_genome}/chr{chr}_{population}_results.tab")
    wildcard_constraints:
        population = '|'.join([population for population in populations if population != 'AA'
                               and population != 'AOUAFR' and population != 'AOUEUR' and population != 'AOUNA'])
    params:
        ibdmix_dir = ibdmix_directory,
        archaic_error_rate = archaic_error_rate,
        max_modern_error_rate = max_modern_error_rate,
        modern_error_proportion = modern_error_proportion,
        min_minor_allele_count = min_minor_allele_count,
        lod_threshold = lod_threshold
    shell:
        "{params.ibdmix_dir}ibdmix -g {input.gt} -o {output} -s <( grep -w -v -f <(cut -f2 {input.kin}) {input.samples}) "
        "-r {input.mask} -m {params.min_minor_allele_count} -a {params.archaic_error_rate} "
        "-e {params.max_modern_error_rate} -c {params.modern_error_proportion} -d {params.lod_threshold}"

use rule run_ibdmix_1kgp as run_ibdmix_aou with:
    input:
        gt = results_path + ibdmix_genotypes_path + "chr{chr}_{archaic_genome}_ibdmix_genotypes_aou",
        samples = data_path + "{population}_sample_ids.txt",
        mask = data_path + mask_path + "chr{chr}.regions_to_exclude.{archaic_genome}.bed",
        kin=related_samples
    output:
        temp(results_path + "ibdmix_{archaic_genome}/chr{chr}_{population}_results.tab")
    wildcard_constraints:
        population = "AOUAFR|AOUEUR|AOUNA"

def get_gt_target(wildcards):
    if dataset == 'test':
        gt = results_path + ibdmix_genotypes_path + f"chr{wildcards.chr}_{wildcards.archaic_genome}_ibdmix_genotypes_target"
    else:
        gt = results_path + ibdmix_genotypes_path + f"chr{wildcards.chr}_{wildcards.archaic_genome}_ibdmix_genotypes_aou"
    return gt

use rule run_ibdmix_1kgp as run_ibdmix_target with:
    input:
        gt=get_gt_target,
        samples = data_path + "AA_sample_ids.txt",
        mask = data_path + mask_path + "chr{chr}.regions_to_exclude.{archaic_genome}.bed",
        kin=related_samples
    output:
        temp(results_path + "ibdmix_{archaic_genome}/chr{chr}_AA_results.tab")

def map_population_to_superpopulation(wildcards):
    return pop_superpop_mapping[wildcards.population]

# format IBDmix results
rule format_ibd_results:
    input:
        results_path + "ibdmix_{archaic_genome}/chr{chr}_{population}_results.tab"
    output:
        results_path + "ibdmix_{archaic_genome}/chr{chr}_{population}_formatted.tab"
    params:
        superpop = map_population_to_superpopulation,
        min_segment_length= minimum_length,
        lod_threshold= lod_threshold
    shell:
        "tail -n+2 {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($4 - $3 >= {params.min_segment_length} && "
        "$5 >= {params.lod_threshold}) print \"chr\"$2,$3,$4,$5,$1,\"{wildcards.population}\","
        "\"{params.superpop}\"}}' > {output}"


def get_input_merge_ibd_results(wildcards):
    return [results_path + f"ibdmix_{wildcards.archaic_genome}/chr{chrom}_{population}_formatted.tab"
            for chrom in chromosomes for population in populations]

# merge IBDmix results for all samples
rule merge_ibd_results:
    input:
        get_input_merge_ibd_results
    output:
        temp(results_path + "ibdmix_{archaic_genome}/" + "ibdmix_results_combined_" +
             str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed")
    shell:
        "cat {input} | sort -k1,1 -k2,2n > {output}"

# mask Denisovan segments found in Africa in Neanderthal call set
rule mask_denisovan_segments_in_neanderthal_set:
    input:
        neanderthal = results_path + "ibdmix_{archaic_genome}/ibdmix_results_combined_" +
                      str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed",
        denisovan = f"{results_path}ibdmix_{denisovan_genome}/ibdmix_results_combined_" +
                    f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
             str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes) # restrict wildcard.archaic_genome to Neanderthal genomes
    params:
        bedtools = bedtools_path,
        min_segment_length = minimum_length,
        lod_threshold = lod_threshold,
        afr_pop = "ESN|GWD|LWK|MSL|YRI"
    shell:
        "{params.bedtools} intersect -v -a {input.neanderthal} "
        "-b <(awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($6~/({params.afr_pop})/) print $0}}' {input.denisovan}) | "
        "awk -F \"\\t\" 'BEGIN{{OFS=\"\\t\"}}{{if ($3 - $2 >= {params.min_segment_length}) print $0}}' > {output}"

rule mask_segments_in_AFR:
    input:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_afr_masked.bed"
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} intersect -a {input} -b <(grep AFR {input}) -v > {output}"

# # mask segments found in AA or AFR but not in EUR or EAS
# rule create_EUR_introgression_reference_panel:
#     input:
#         results_path + "tmp_ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#         str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed"
#     output:
#         temp(results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#              str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_EUR_reference.bed")
#     params:
#         bedtools=bedtools_path
#     shell:
#         "grep -E 'EUR' {input} | bedtools merge > {output}"
#
# rule filter_introgressed_segments_by_EUR_reference_panel:
#     input:
#         target=results_path + "tmp_ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#         str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed",
#         panel=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#         str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_EUR_reference.bed"
#     output:
#         results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#         str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed"
#     params:
#         bedtools=bedtools_path,
#         filter_pops = 'AMR|AFR',
#         non_filter_pops = 'EAS|EUR'
#     shell:
#         "grep -E \"{params.filter_pops}\" {input.target} |"
#         " {params.bedtools} intersect -a - -b {input.panel} -f 1 -u > {output}.tmp; "
#         "grep -E \"{params.non_filter_pops}\" {input.target} >> {output}.tmp; "
#         "sort -k1,1 -k2,2n {output}.tmp > {output}; rm {output}.tmp"
#
# rule remove_outlier_admixed_individuals:
#     input:
#         results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#         str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed"
#     output:
#         filtered = results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#                    str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_removed_admixed_outlier_individuals.bed",
#         outliers = results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
#                    str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_admixed_outlier_individuals.bed",
#         outlier_ids = results_path + "ibdmix_{archaic_genome}/admixed_outlier_iids.txt"
#     run:
#         ibdmix = pd.read_csv(input[0], sep='\t', names=['chrom', 'start', 'end', "LOD", "IID",
#                                                         "pop", "super_pop"])
#         ibdmix['length'] = (ibdmix.end - ibdmix.start)
#         amounts = ibdmix[ibdmix['pop'] == 'AA'].groupby('IID').sum().loc[:, ['length']] / 1e6
#         # # remove outlier individuals
#         amounts = amounts[(amounts['length'] <= np.percentile(amounts['length'].values, 2.5)) |
#                           (amounts['length'] >= np.percentile(amounts['length'].values, 97.5))]
#         ibdmix.drop('length', axis=1, inplace=True)
#         ibdmix_filtered = ibdmix[~np.isin(ibdmix.IID, amounts.index.values)]
#         ibdmix_outliers = ibdmix[np.isin(ibdmix.IID, amounts.index.values)]
#         ibdmix_filtered.to_csv(output.filtered, sep='\t', index=False, header=False)
#         ibdmix_outliers.to_csv(output.outliers, sep='\t', index=False, header=False)
#         with open(output.outlier_ids, 'w') as out:
#             for iid in amounts.index.values:
#                 out.write(f'{iid}\n')
#         out.close()



rule calculate_introgressed_coverage_per_window_and_superpopulation:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_afr_masked.bed",
        genomefile=data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes) # restrict wildcard.archaic_genome to Neanderthal genomes
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} intersect -a {input.genomefile} -b {input.ibdmix} -wo | sort -k1,1 -k2,2n -k3,3n -k10 | "
        "{params.bedtools} groupby -g 1,2,3,10 -c 11 -o sum > {output}"

rule calculate_introgressed_coverage_per_window_and_individual:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_afr_masked.bed",
        genomefile=data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) +
        "LOD_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes), # restrict wildcard.archaic_genome to Neanderthal genomes
        windowsize=windowsize,
        stepsize=stepsize
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} intersect -a {input.genomefile} -b {input.ibdmix} -wo | sort -k1,1 -k2,2n -k3,3n -k8 | "
        "{params.bedtools} groupby -g 1,2,3,8 -c 11 -o sum > {output}"

rule split_introgressed_coverage_per_window_and_superpopulation:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_afr_masked_coverage_per_window{windowsize}_s_{stepsize}.bed",
    output:
        results_path + "ibdmix_{archaic_genome}/{superpopulation}_ibdmix_results_masked_denisovan_combined_" +
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_afr_masked_frequencies_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes) # restrict wildcard.archaic_genome to Neanderthal genomes
    shell:
        "grep -w {wildcards.superpopulation} {input.ibdmix} > {output}"

use rule calculate_introgressed_coverage_per_window_and_superpopulation as calculate_introgressed_coverage_per_window_and_superpopulation_all_segments with:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +\
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed",
        genomefile=data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" + \
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_coverage_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes) # restrict wildcard.archaic_genome to Neanderthal genomes
    params:
        bedtools = bedtools_path

use rule calculate_introgressed_coverage_per_window_and_individual as calculate_introgressed_coverage_per_window_and_individual_all_segments with:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +\
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed",
        genomefile=data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +\
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) +\
        "LOD_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes), # restrict wildcard.archaic_genome to Neanderthal genomes
        windowsize=windowsize,
        stepsize=stepsize
    params:
        bedtools = bedtools_path

use rule split_introgressed_coverage_per_window_and_superpopulation as split_introgressed_coverage_per_window_and_superpopulation_all_segments with:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +\
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_coverage_per_window{windowsize}_s_{stepsize}.bed",
    output:
        results_path + "ibdmix_{archaic_genome}/{superpopulation}_ibdmix_results_masked_denisovan_combined_" +\
        str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD_frequencies_per_window{windowsize}_s_{stepsize}.bed"
    wildcard_constraints:
        archaic_genome="|".join(neanderthal_genomes) # restrict wildcard.archaic_genome to Neanderthal genomes

def get_masks(wildcards):
    return [data_path + mask_path + f"chr{chrom}.regions_to_exclude.{wildcards.archaic_genome}.bed"
            for chrom in chromosomes]

def get_mask_pattern(wildcards):
    return data_path + mask_path + 'chr{chrom}.regions_to_exclude.' + f'{wildcards.archaic_genome}.bed'

rule find_putatively_selected_neanderthal_segments_african_americans:
    input:
        masks = get_masks,
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_window{windowsize}_s_{stepsize}"
             + ".{archaic_genome}.bed" for chrom in chromosomes],
        nea_amr=results_path + "ibdmix_{archaic_genome}/AMR_ibdmix_results_masked_denisovan_combined_" +
                f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD_afr_masked_frequencies_per_window{windowsize}_s_{stepsize}.bed",
        nea_eur=results_path + "ibdmix_{archaic_genome}/EUR_ibdmix_results_masked_denisovan_combined_" +
                f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD_afr_masked_frequencies_per_window{windowsize}_s_{stepsize}.bed",
        nea_eas=results_path + "ibdmix_{archaic_genome}/EAS_ibdmix_results_masked_denisovan_combined_" +
                f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD_afr_masked_frequencies_per_window{windowsize}_s_{stepsize}.bed",
        ibdmix_all=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
                   f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD_afr_masked.bed",
        ibdmix_all_coverage = results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
                              str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) +
                              f"LOD_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed",
        genomefile=data_path + reference_path + 'genomefile_hg38.bed',
        windowed_genomefile=data_path + reference_path + f'hg38_windowed_w_{windowsize}_s_{stepsize}.bed',

        genetic_maps = expand(data_path + genetic_map_path + 'genetic_map_Hg38_chr{chrom}.txt', chrom=chromosomes),
        chrom_sizes=data_path + reference_path + 'hg38.chrom.sizes.bed'
    output:
        selected=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments.bed",
        not_selected=[results_path + "ibdmix_{archaic_genome}/" +
                      f"AMR_putatively_not_selected_control_neanderthal_segments_{n}.bed"
                      for n in range(bootstrap_reps)],
        pvals=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
                              str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) +
                              f"LOD_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_pvalues.bed",
        exp=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
              str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) +
              f"LOD_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_expectations.bed",
    params:
        min_length = minimum_length_selected,
        max_eas_anc = 1 - min_afr_eur_component_combined,
        window_size = windowsize,
        step_size=stepsize,
        alpha = alpha,
        min_expectation = min_expectation,
        not_selected_prefix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_not_selected_control_neanderthal_segments",
        reps = bootstrap_reps,
        mask_pattern = get_mask_pattern
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes)  # restrict wildcard.archaic_genome to Neanderthal genomes
    threads: 64
    shell:
        "scripts/identify_putatively_selected_introgressed_segments_in_african_americans.py -l {input.lai} "
        "--nea_amr {input.nea_amr} --nea_eur {input.nea_eur} --nea_eas {input.nea_eas} "
        "--ibdmix_all {input.ibdmix_all} --ibdmix_all_coverage {input.ibdmix_all_coverage} "
        "--min_length {params.min_length} --max_eas_anc {params.max_eas_anc} --min_expectation {params.min_expectation} "
        "-g {input.genomefile} -os {output.selected} -s {params.step_size} "
        "-ons {params.not_selected_prefix} -oip {output.pvals} -oie {output.exp} -w {params.window_size} "
        "--alpha {params.alpha} --reps {params.reps} --threads {threads} --genetic_maps {input.genetic_maps} "
        "--chrom_sizes {input.chrom_sizes} --windowed_genome_file {input.windowed_genomefile} "
        "--mask_pattern {params.mask_pattern}"

def get_lai_pattern(wildcards):
    return results_path + flare_output + "african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}" + \
           f".{wildcards.archaic_genome}.bed"

rule find_neanderthal_introgression_deserts_AMR:
    input:
        masks=get_masks,
        coverage= results_path + 'ibdmix_{archaic_genome}/AMR_ibdmix_results_masked_denisovan_combined_50kb_4.0' +
                 f'LOD_frequencies_per_window{windowsize}_s_{stepsize}.bed',
        coverage_masked = results_path + 'ibdmix_{archaic_genome}/AMR_ibdmix_results_masked_denisovan_combined_50kb_4.0' +
                          f'LOD_afr_masked_frequencies_per_window{windowsize}_s_{stepsize}.bed',
        genomefile = data_path + reference_path + "genomefile_hg38.bed",
        ibdmix_all=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
                   f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD.bed",
        ibdmix_afr_masked = results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
                   f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD_afr_masked.bed",
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}" +
             ".{archaic_genome}.bed" for chrom in chromosomes for phase in [0, 1]],
        # cw_afr = results_path+ "ibdmix_{archaic_genome}/AFR_introgression_frequencies_and_rank_callable_windows.bed",
        cw_eur= results_path+ "ibdmix_{archaic_genome}/EUR_introgression_frequencies_and_rank_callable_windows_afr_masked.bed",
        cw_eas= results_path+ "ibdmix_{archaic_genome}/EAS_introgression_frequencies_and_rank_callable_windows_afr_masked.bed"
    output:
        deserts=results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts.bed",
        controls=[results_path + "ibdmix_{archaic_genome}/" + "AMR_introgression_deserts_control_segments_"
                  + f"{n}.bed" for n in range(bootstrap_reps)],
        controls_new=[results_path + "ibdmix_{archaic_genome}/" + "AMR_introgression_deserts_new_control_segments_"
                      + f"{n}.bed" for n in range(bootstrap_reps)],
        ranks=results_path + "ibdmix_{archaic_genome}/AMR_introgression_frequencies_and_rank_callable_windows.bed",
        ranks_masked=results_path + "ibdmix_{archaic_genome}/AMR_introgression_frequencies_and_rank_callable_windows_afr_masked.bed",
        deserts_new=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts.bed",
        deserts_new_pvalues=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_pvalues.bed"
    params:
        mask_pattern = get_mask_pattern,
        max_masked = max_masked,
        stride = stride,
        window_sizes = ' '.join([str(w) for w in window_sizes]),
        reps= bootstrap_reps,
        control_prefix = results_path+ "ibdmix_{archaic_genome}/AMR_introgression_deserts_control_segments",
        control_new_prefix= results_path+ "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments",
        super_pop = 'AMR',
        lai_pattern = get_lai_pattern,
        tmp_dir = tmp_dir
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes)  # restrict wildcard.archaic_genome to Neanderthal genomes
    threads: 64
    resources:
        mem_mb=128 * 1000,
        load=5
    shell:
        'scripts/find_introgression_deserts.py --introgression_coverage {input.coverage} '
        '--introgression_coverage_masked {input.coverage_masked} --stride {params.stride} '
        '--max_masked {params.max_masked} --window_sizes {params.window_sizes} --reps {params.reps} '
        '--genomefile {input.genomefile} --mask_pattern {params.mask_pattern} -o {output.deserts} '
        '--ibdmix {input.ibdmix_all} --ibdmix_masked {input.ibdmix_afr_masked} -or {params.control_prefix} '
        '--threads {threads} --output_ranks {output.ranks} --output_ranks_masked {output.ranks_masked} '
        '--super_pop {params.super_pop} --lai {params.lai_pattern} --callable_windows_references {input.cw_eur} '
        '{input.cw_eas} --callable_windows_references_labels EUR EAS '
        '--unique_deserts {output.deserts_new} --unique_deserts_pvalues {output.deserts_new_pvalues} '
        '-orn {params.control_new_prefix} --tmp_dir {params.tmp_dir}'

rule find_neanderthal_introgression_deserts:
    input:
        masks=get_masks,
        coverage= results_path + 'ibdmix_{archaic_genome}/{superpopulation}_ibdmix_results_masked_denisovan_combined_50kb_4.0' +
                 f'LOD_frequencies_per_window{windowsize}_s_{stepsize}.bed',
        coverage_masked = results_path + 'ibdmix_{archaic_genome}/{superpopulation}_ibdmix_results_masked_denisovan_combined_50kb_4.0' +
                          f'LOD_afr_masked_frequencies_per_window{windowsize}_s_{stepsize}.bed',
        genomefile = data_path + reference_path + "genomefile_hg38.bed",
        ibdmix_all=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
                   f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD.bed",
        ibdmix_afr_masked= results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
        f"{int(minimum_length / 1000)}kb_{lod_threshold}LOD_afr_masked.bed",
    output:
        deserts=results_path + "ibdmix_{archaic_genome}/{superpopulation}_introgression_deserts.bed",
        ranks=results_path+ "ibdmix_{archaic_genome}/{superpopulation}_introgression_frequencies_and_rank_callable_windows.bed",
        ranks_masked=results_path+ "ibdmix_{archaic_genome}/{superpopulation}_introgression_frequencies_and_rank_callable_windows_afr_masked.bed"
    params:
        mask_pattern = get_mask_pattern,
        max_masked = max_masked,
        stride = stride,
        window_sizes = ' '.join([str(w) for w in window_sizes]),
        super_pop = '{superpopulation}',
        tmp_dir= tmp_dir
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes), # restrict wildcard.archaic_genome to Neanderthal genomes,
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
        '--ibdmix {input.ibdmix_all} --ibdmix_masked {input.ibdmix_afr_masked} --tmp_dir {params.tmp_dir}'

# # RUN ARCHAIC HMM FROM SKOV ET AL.
# rule convert_vcf_to_bcf:
#     input:
#         data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.vcf.gz'
#     output:
#         data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf'
#     params:
#         bcftools=bcftools_path
#     shell:
#         "{params.bcftools} view -l1 -Ob -o {output} {input}"
#
# rule index_bcf:
#     input:
#         data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf'
#     output:
#         data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf.csi'
#     params:
#         bcftools=bcftools_path
#     shell:
#         "{params.bcftools} index {input}"
#
# def get_ingroup_sample_ids(wildcards):
#     non_afr_populations = [key for key, val in pop_superpop_mapping.items() if val != "AFR"]
#     return " ".join([data_path + f"{population}_sample_ids.txt" for population in non_afr_populations])
#
# def get_outgroup_sample_ids(wildcards):
#     afr_populations = [key for key, val in pop_superpop_mapping.items() if val == "AFR"]
#     return " ".join([data_path + f"{population}_sample_ids.txt" for population in afr_populations])
#
# rule create_individual_json:
#     input:
#         sample_ids=expand(data_path + "{population}_sample_ids.txt", population=populations),
#         kin=related_samples
#     output:
#         data_path + "individuals_hmmix.json"
#     params:
#         ingroup_ids=get_ingroup_sample_ids,
#         outgroup_ids=get_outgroup_sample_ids
#     shell:
#         "scripts/create_indviduals_json.py --ingroup {params.ingroup_ids} --outgroup {params.outgroup_ids} "
#         "--out {output} --kin {input.kin}"
#
# rule get_sites_to_include:
#     input:
#         mask=data_path + mask_path + "chr{chr}.regions_to_exclude.Vindija33.19.bed",
#         chrom_sizes=data_path + reference_path + hg38_chrom_sizes_url.split('/')[-1] + '.bed'
#     output:
#         data_path + mask_path + "chr{chr}.regions_to_include.Vindija33.19.bed"
#     params:
#         bedtools = bedtools_path,
#         chrom = "{chr}"
#     shell:
#         "{params.bedtools} subtract -a <(grep -w chr{params.chrom} {input.chrom_sizes} | sed 's/chr//' ) "
#         "-b {input.mask} | sed 's/^/chr/' > {output}"
#
# rule merge_mask_hmmix:
#     input:
#         mask = expand(data_path + mask_path + "chr{chr}.regions_to_include.Vindija33.19.bed",chr=chromosomes),
#     output:
#         data_path + mask_path + "hmmix.regions_to_include.Vindija33.19.bed"
#     shell:
#         "cat {input} | sort -k1,1 -k2,2n > {output}"
#
# rule split_reference_genome_by_chrom:
#     input:
#         data_path + reference_path + "hg38.fa"
#     output:
#         expand(data_path + reference_path + "chr{chrom}_hg38.fa", chrom=chromosomes)
#     params:
#         prefix = data_path + reference_path,
#         suffix = "_hg38.fa"
#     shell:
#         "scripts/split_reference_fasta_by_chromosome.py --input {input} --prefix {params.prefix} "
#         "--suffix {params.suffix}"
#
# def get_bcf(wildcards):
#     return ','.join([data_path + merged_datasets_path  + f'reference_target_individuals_chr{chrom}.bcf' for chrom in chromosomes])
#
# def get_bcf_index(wildcards):
#     return ','.join([data_path + merged_datasets_path  + f'reference_target_individuals_chr{chrom}.bcf.csi' for chrom in chromosomes])
#
# def get_ancestral_sequences(wildcards):
#     return ','.join([data_path + human_ancestral_sequence_path + f'homo_sapiens_ancestor_{chrom}.fa' for chrom in chromosomes])
#
# def get_reference_sequences(wildcards):
#     return ",".join([data_path + reference_path + f"chr{chrom}_hg38.fa" for chrom in chromosomes])
#
# rule create_outgroup:
#     input:
#         bcf=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf', chr=chromosomes),
#         csi=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf.csi', chr=chromosomes),
#         individuals=data_path + "individuals_hmmix.json",
#         mask=data_path + mask_path + "hmmix.regions_to_include.Vindija33.19.bed",
#         ancestral=expand(data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_{chr}.fa', chr=chromosomes),
#         reference=expand(data_path + reference_path + "chr{chrom}_hg38.fa", chrom=chromosomes)
#     output:
#         data_path + "outgroup_hmmix.txt"
#     conda:
#         "../envs/archaic_hmm.yaml"
#     params:
#         bcf=get_bcf,
#         ancestral=get_ancestral_sequences,
#         reference=get_reference_sequences
#     shell:
#         "hmmix create_outgroup -ind={input.individuals} -vcf={params.bcf} -out={output} -weights={input.mask} "
#         "-refgenome={params.reference} -ancestral={params.ancestral}"
#
# rule estimate_mutation_rate:
#     input:
#         outgroup=data_path + "outgroup_hmmix.txt",
#         mask=data_path + mask_path + "hmmix.regions_to_include.Vindija33.19.bed",
#     output:
#         data_path + "mutationrate_hmmix.bed"
#     conda:
#         "../envs/archaic_hmm.yaml"
#     shell:
#         "hmmix mutation_rate -outgroup={input.outgroup} -window_size=1000000 -out={output} -weights={input.mask}"
#
#
# checkpoint create_ingroup:
#     input:
#         bcf=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf', chr=chromosomes),
#         csi=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bcf.csi', chr=chromosomes),
#         individuals=data_path + "individuals_hmmix.json",
#         outgroup=data_path + "outgroup_hmmix.txt",
#         mask=data_path + mask_path + "hmmix.regions_to_include.Vindija33.19.bed",
#         ancestral=expand(data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_{chr}.fa',chr=chromosomes),
#     output:
#         directory(results_path + "hmmix/")
#     params:
#         bcf=get_bcf,
#         ancestral=get_ancestral_sequences,
#         output_prefix=results_path + "hmmix/obs"
#     conda:
#         "../envs/archaic_hmm.yaml"
#     shell:
#         "hmmix create_ingroup  -ind={input.individuals} -vcf={params.bcf} -out={params.output_prefix} "
#         "-outgroup={input.outgroup} -weights={input.mask} -ancestral={params.ancestral}"
#
# rule train_hmm:
#     input:
#         obs=results_path + "hmmix/obs.{sample}.txt",
#         mutrates=data_path + "mutationrate_hmmix.bed",
#         mask=data_path + mask_path + "hmmix.regions_to_include.Vindija33.19.bed"
#     output:
#         results_path + "hmmix/trained.{sample}.txt"
#     conda:
#         "../envs/archaic_hmm.yaml"
#     shell:
#         "hmmix train  -obs={input.obs} -mutrates={input.mutrates} -out={output} -haploid -weights={input.mask}"
#
# rule decode_hmm:
#     input:
#         obs=results_path + "hmmix/obs.{sample}.txt",
#         mutrates=data_path + "mutationrate_hmmix.bed",
#         trained=results_path + "hmmix/trained.{sample}.txt",
#         mask=data_path + mask_path + "hmmix.regions_to_include.Vindija33.19.bed"
#     output:
#         multiext(results_path + "hmmix/{sample}.decoded", ".hap1.txt", ".hap2.txt")
#     params:
#         out_prefix = results_path + "hmmix/{sample}.decoded"
#     conda:
#         "../envs/archaic_hmm.yaml"
#     shell:
#         "hmmix decode -obs={input.obs} -mutrates={input.mutrates} -param={input.trained} -haploid "
#         "-out={params.out_prefix} -weights={input.mask}"
#
# rule format_hmm_output:
#     input:
#         results_path + "hmmix/{sample}.decoded.hap{phase}.txt"
#     output:
#         results_path + "hmmix/{sample}.decoded.formatted.hap{phase}.txt"
#     params:
#         sample="{sample}",
#         phase="{phase}"
#     shell:
#         "grep Archaic {input} | "
#         "awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($6 >= 0.8) print $1, $2, $3, $6, \"{params.sample}\", \"{params.phase}\"}}' "
#         "> {output}"
#
# rule merge_phases:
#     input:
#         phase1=results_path + "hmmix/{sample}.decoded.formatted.hap1.txt",
#         phase2=results_path + "hmmix/{sample}.decoded.formatted.hap2.txt"
#     output:
#         results_path + "hmmix/{sample}.decoded.formatted.txt"
#     shell:
#         "cat {input} | sort -k1,1 -k2,2n > {output}"
#
# def get_all_samples(wildcards):
#     checkpoints.create_ingroup.get(**wildcards)
#     samples = glob_wildcards(results_path + "hmmix/obs.{sample}.txt").sample
#     return expand(results_path + "hmmix/{sample}.decoded.formatted.txt", sample=samples)
#
# rule aggregate_samples_hmmix:
#     input:
#         get_all_samples
#     output:
#         results_path + "introgressed_segments_hmmix.bed"
#     shell:
#         "cat {input} | sort -k1,1 -k2,2n > {output}"
