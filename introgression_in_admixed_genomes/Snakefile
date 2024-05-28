import numpy as np
from snakemake.remote.GS import RemoteProvider as GSRemoteProvider
np.random.seed(42)

GS = GSRemoteProvider(project="terra-vpc-sc-373c109f")
GS_PREFIX = "fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/acaf_threshold_v7.1/vcf"

configfile: "config/config.yaml"

related_samples = config['related_samples']
all_of_us_vcf = config['all_of_us_vcf']
phase3_1KG_base_url = config['phase3_1KG_base_url']
hgdp_wgs_url = config['hgdp_wgs_url']
hgdp_mayan_samples = config['hgdp_mayan_samples']
hgdp_pima_samples = config['hgdp_pima_samples']
beagle_url = config['beagle_url']
accessible_genome_mask_url = config['accessible_genome_mask_url']
segmental_duplications_url = config['segmental_duplications_url']
reference_gap_url = config['reference_gap_url']
knownGene_url = config['knownGene_url']
knownToEnsembl_url = config['knownToEnsembl_url']
gerp_url = config['gerp_url']
phastCons_url = config['phastCons_url']
phyloP_url = config['phyloP_url']
recomb_rate_url = config['recomb_rate_url']
bstat_url = config['bstat_url']
gnomad_url = config['gnomad_url']
clinvar_url = config['clinvar_url']
anno_fields = config['anno_fields']
encode_annotation_url = config['encode_annotation_url']
simple_repeat_map_url = config['simple_repeat_map_url']
hg38_chrom_sizes_url = config['hg38_chrom_sizes_url']
hg38_fasta_url = config['hg38_fasta_url']
hg38_cytobands_url = config['hg38_cytobands_url']
string_url = config['string_url']
gwas_efo_trait_mapping_url = config['gwas_efo_trait_mapping_url']
gwas_catalog_url = config['gwas_catalog_url']
dbSNP_url = config['dbSNP_url']
whole_genome_alignments_url = config['whole_genome_alignments_url']
seqbility_url = config['seqbility_url']
genetic_map_plink_url = config['hapmap_genetic_map_plink_url']
genetic_map_url = config['hapmap_genetic_map_url']
ensembl_uniprot_url = config['ensembl_uniprot_url']
gtex_url = config['gtex_url']
species = config['species']
archaic_genomes_mpg_base_url = config['archaic_genomes_mpg_base_url']
archaic_genomes_mpg_base_url_altai_old = config['archaic_genomes_mpg_base_url_altai_old']
archaic_genomes_mpg_url_filter_bed_altai_old = config["archaic_genomes_mpg_url_filter_bed_altai_old"]
hg19_to_hg38_chain_url = config['hg19_to_hg38_chain_url']
hg38_to_hg19_chain_url = config['hg38_to_hg19_chain_url']
url_boosting_scores_eur_complete = config['url_boosting_scores_eur_complete']
url_boosting_scores_eur_recent_complete = config['url_boosting_scores_eur_recent_complete']
url_boosting_scores_eur_ancient_complete = config['url_boosting_scores_eur_ancient_complete']
url_boosting_scores_afr_complete = config['url_boosting_scores_afr_complete']
url_boosting_scores_afr_recent_complete = config['url_boosting_scores_afr_recent_complete']
url_boosting_scores_afr_ancient_complete = config['url_boosting_scores_afr_ancient_complete']
url_boosting_scores_eas_complete = config['url_boosting_scores_eas_complete']
url_boosting_scores_eas_recent_complete = config['url_boosting_scores_eas_recent_complete']
url_boosting_scores_eas_ancient_complete = config['url_boosting_scores_eas_ancient_complete']

url_boosting_scores_eur_incomplete = config['url_boosting_scores_eur_incomplete']
url_boosting_scores_eur_recent_incomplete = config['url_boosting_scores_eur_recent_incomplete']
url_boosting_scores_eur_ancient_incomplete = config['url_boosting_scores_eur_ancient_incomplete']
url_boosting_scores_afr_incomplete = config['url_boosting_scores_afr_incomplete']
url_boosting_scores_afr_recent_incomplete = config['url_boosting_scores_afr_recent_incomplete']
url_boosting_scores_afr_ancient_incomplete = config['url_boosting_scores_afr_ancient_incomplete']
url_boosting_scores_eas_incomplete = config['url_boosting_scores_eas_incomplete']
url_boosting_scores_eas_recent_incomplete = config['url_boosting_scores_eas_recent_incomplete']
url_boosting_scores_eas_ancient_incomplete = config['url_boosting_scores_eas_ancient_incomplete']

archaic_genomes_mpg_base_url_chagyrskaya = config['archaic_genomes_mpg_base_url_chagyrskaya']
archaic_genomes = config['archaic_genomes']
neanderthal_genomes = config['neanderthal_genomes']
denisovan_genome = config['denisovan_genome']
populations = config['populations']
pop_superpop_mapping = config['pop_superpop_mapping']
ibdmix_directory = config['ibdmix_directory']
simulations_path = config['simulations_path_base']
bcftools_path = config["bcftools_path"]
bedtools_path = config["bedtools_path"]
bwa_path = config['bwa_path']
flare_path = config['flare_path']
ensembl_path = config['ensembl_path']
rye_path = config['rye_path']
data_path = config['data_path']
ibdmix_genotypes_path = config['ibdmix_genotypes_path']
results_path = config['results_path']
ldsc_path = config['ldsc_path']
tmp_dir = config['tmp_dir']
ldsc_plink_url = config['ldsc_plink_url']
ldsc_freq_url = config['ldsc_freq_url']
ldsc_weights_url = config['ldsc_weights_url']
ldsc_baseline_model_url = config['ldsc_baseline_model_url']
gwas_summary_stats_url = config['gwas_summary_stats_url']
mask_path = config['mask_path']
whole_genome_alignments_path = config['whole_genome_alignments_path']
reference_path = config['reference_path']
seqbility_tmp = config['seqbility_tmp']
merged_datasets_path = config['merged_datasets_path']
flare_output = config['flare_output_path']
genetic_map_path = config["genetic_map_path"]
human_ancestral_sequence_path = config['human_ancestral_sequence_path']
maf = config['maf']
geno = config['geno']
minimum_length = config['minimum_length']
minimum_length_selected = config['minimum_length_selected']
lod_threshold = config['lod_threshold']
archaic_error_rate = config['archaic_error_rate']
max_modern_error_rate = config['max_modern_error_rate']
modern_error_proportion = config['modern_error_proportion']
min_minor_allele_count = config['min_minor_allele_count']
K = config['K']
M = config['M']
min_afr_component = config['min_afr_component']
max_afr_component = config['max_afr_component']
min_eur_component = config['min_eur_component']
max_eur_component = config['max_eur_component']
n_replicates = np.arange(config['simulation_replicates'])
min_afr_eur_component_combined = config['min_afr_eur_component_combined']
max_masked = config['max_masked']
stride = config['stride']
window_sizes = config['window_sizes']
chromosomes = np.arange(1, 23).tolist()
windowsize = config['windowsize']
stepsize = config['stepsize']
alpha = config['alpha']
min_expectation = config['min_expectation']
generation_time = config['generation_time']
mutation_rate = config['mutation_rate']
arg_dir = config['arg_dir']
effective_population_size = config['effective_population_size']
n_sites_idat = config['n_sites_idat']
windowsize_idat = config['windowsize_idat']
stepsize_idat = config['stepsize_idat']
dist_step_size = config['dist_step_size']
min_cov_dat = config['min_cov_dat']
mem_gb_cpu = config['mem_gb_cpu']
bootstrap_reps = config['bootstrap_reps']
dataset = config['dataset']
ldlink_token = config['ldlink_token']

include: "snakefiles/download_and_prepare_data.smk"
include: "snakefiles/generate_mappability_mask.smk"
include: "snakefiles/determine_ancestry.smk"
include: "snakefiles/infer_human_ancestral_sequence.smk"
include: "snakefiles/ibdmix.smk"
include: "snakefiles/functional_annotation.smk"
include: "snakefiles/simulations.smk"
include: "snakefiles/estimate_allele_ages.smk"
include: "snakefiles/ldsc.smk"
if dataset == 'test':
    include: "snakefiles/get_target_data_from_1kgp.smk"
    populations = [pop for pop in populations if pop != 'AOUAFR' and pop != 'AOUEUR' and pop != 'AOUNA']
elif dataset == 'AOU':
    include: "snakefiles/phase_all_of_us_data.smk"

rule all:
    input:
        ## generate only input required to run IBDmix
        ## takes 32GB per CPU
        expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.updated_ids{ext}",
              ext=[".bed", ".bim", ".fam"], chr=chromosomes),
        expand(data_path+ "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz",
             archaic_genome=archaic_genomes, chr=chromosomes),
        expand(data_path+ mask_path + "chr{chr}.regions_to_exclude.{archaic_genome}.bed",
             archaic_genome=archaic_genomes, chr=chromosomes),
        expand(data_path+ "{population}_sample_ids.txt", population=[pop for pop  in populations if pop != 'AA']),
        multiext(data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz', ''),
                 '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path + reference_path + 'knownCDS', '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path + reference_path + gerp_url.split('/')[-1].replace('.bw', ''), '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz', ''),
                 '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz', ''),
                 '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', ''),
                 '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38'), '.bed.gz', '.bed.gz.tbi'),
        multiext(data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz', ''),
                '.bed.gz', '.bed.gz.tbi'),
        expand(data_path+ "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100{ext}",
               archaic_genome=archaic_genomes, chr=chromosomes, ext=['.bed.gz', '.bed.gz.tbi']),
        expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated.vcf.gz",
              chr=chromosomes),
        data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz', ''),
        expand(data_path + genetic_map_path + "plink.chr{chr}.GRCh38.map", chr=chromosomes),
        data_path + reference_path + 'genomefile_hg38.bed',
        expand(data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed",
               windowsize=windowsize, stepsize=stepsize),
        expand(data_path + "{population}_sample_ids.txt", population=populations[:-2]),
        expand(data_path + "{superpopulation}_sample_ids.txt", superpopulation=['AFR', 'EUR', 'EAS']),
        data_path + reference_path + ensembl_uniprot_url.split('/')[-1].replace('.gz',''),
        data_path + reference_path + string_url.split('/')[-1].replace('.gz',''),
        expand(data_path + "1000G_phase3/ACB_ASW.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes{ext}",
               ext=[".vcf.gz", ".vcf.gz.tbi"], chr=chromosomes),
        # ## generate IBDmix output that needs to be run on GCP
        # ## takes 8GB per CPU
        expand(results_path + "ibdmix_{archaic_genome}/background_list_of_genes_with_introgression_converted.txt",
               archaic_genome=neanderthal_genomes),
        expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
               archaic_genome=neanderthal_genomes),
        expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
            archaic_genome=neanderthal_genomes, segment_type=[ 'not_selected_control'], n=np.arange(0, bootstrap_reps)),
        expand(data_path + "{superpopulation}_chr{chr}.afreq", superpopulation=['AA', 'AFR', 'EUR', "EAS"],
               chr=chromosomes),
        expand(results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_window{windowsize}_s_{stepsize}.{archaic_genome}.bed",
               chr=chromosomes, windowsize=windowsize, stepsize=stepsize, archaic_genome=neanderthal_genomes),
        expand(results_path + 'ibdmix_{archaic_genome}/iDAT_scores_{chr}.bed', chr=chromosomes, archaic_genome=neanderthal_genomes),
        expand(results_path + 'ibdmix_{archaic_genome}/standardized_iDAT_scores.bed', archaic_genome=neanderthal_genomes),
        expand(results_path +
               "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_pvalues.bed",
               archaic_genome=neanderthal_genomes, windowsize=windowsize, stepsize=stepsize),
        expand(results_path +
               "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_afr_masked_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_expectations.bed",
               archaic_genome=neanderthal_genomes, windowsize=windowsize, stepsize=stepsize),
        expand(results_path + "ibdmix_{archaic_genome}/{super_population}_introgression_frequencies_and_rank_callable_windows_afr_masked.bed",
               archaic_genome=neanderthal_genomes,super_population=['AMR', 'EUR', "EAS"]),
        expand(results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
               archaic_genome=neanderthal_genomes),
        expand(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
               archaic_genome=neanderthal_genomes, n=np.arange(0, bootstrap_reps)),
        #results_path + "nucleotide_diversity_per_10kb.bed.gz",
        # # also run Skovs HMMIX
        # results_path + "introgressed_segments_hmmix.bed",

        # ## post-processing of IBDmix output that can be run on PACE
        # ## need to download:
        # ## results_path + "ibdmix_{archaic_genome}/background_list_of_genes_with_introgression_converted.txt"
        # ## results_path + 'ibdmix_{archaic_genome}/iDAT_scores_{chr}.bed'
        # ## results_path + 'ibdmix_{archaic_genome}/standardized_iDAT_scores.bed'
        # ## data_path + "{superls d}_chr{chr}.afreq"
        # ## results_path+ flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase0.{archaic_genome}.bed"
        # ## results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase1.{archaic_genome}.bed"
        # ## results_path + ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_coverage_per_individual_and_per_window10000.bed
        # ## results_path + ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_coverage_per_individual_and_per_window10000_pvalues.bed
        # ## results_path + ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_50kb_4.0LOD_coverage_per_window10000.bed
        # ## results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed"
        # ## results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed"
        # ## expand(results_path + "ibdmix_{archaic_genome}/{super_population}_introgression_frequencies_and_rank_callable_windows.bed",
        # ##       archaic_genome=neanderthal_genomes,super_population=['AFR', 'AMR', 'EUR', "EAS"]),
        # ## expand(results_path + "ibdmix_{archaic_genome}/{superpopulation}_novel_introgression_deserts_iDAT_annotated.bed",
        # ##       archaic_genome=neanderthal_genomes, superpopulation=['AFR', 'AMR', 'EUR', "EAS"]),
        # ## expand(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        # ##      archaic_genome=neanderthal_genomes, n=np.arange(0, bootstrap_reps)),
        # ## results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_pvalues.bed"
        # ## results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts.bed"
        # ## results_path+ "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed"
        # ## results_path + "nucleotide_diversity_per_10kb.bed.gz"
        # ## takes 32GB per CPU
        # data_path+ reference_path + 'gwas_catalog_trait_mapping.tab',
        # expand(results_path + "ibdmix_{archaic_genome}/foreground_list_of_genes_with_{selection}_selected_introgression_converted.txt",
        #        archaic_genome=neanderthal_genomes, selection=['pos', 'neg']),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_independent_gwas_hits.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_independent_gwas_hits.bed",
        #        archaic_genome=neanderthal_genomes, segment_type=['not_selected_control'],
        #        n=np.arange(0, bootstrap_reps)),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_independent_gwas_hits_parent_terms.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_independent_gwas_hits_parent_terms.bed",
        #        archaic_genome=neanderthal_genomes, segment_type=['not_selected_control'],
        #        n=np.arange(0, bootstrap_reps)),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_eQTLs.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_eQTLs.bed",
        #        archaic_genome=neanderthal_genomes, segment_type=['not_selected_control'],
        #        n=np.arange(0, bootstrap_reps)),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_annotated.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_annotated.bed",
        #        segment_type=['not_selected_control'], n=np.arange(0,bootstrap_reps), archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_annotated.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_foreground_list_of_genes_in_deserts_converted.txt",
        #        archaic_genome=neanderthal_genomes,),
        # expand(results_path+ "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_independent_gwas_hits.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path+ "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_independent_gwas_hits_parent_terms.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_eQTLs.bed",
        #        archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_annotated.bed",
        #        n=np.arange(0, bootstrap_reps), archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_independent_gwas_hits.bed",
        #        n=np.arange(0, bootstrap_reps), archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_independent_gwas_hits_parent_terms.bed",
        #        n=np.arange(0, bootstrap_reps), archaic_genome=neanderthal_genomes),
        # expand(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_eQTLs.bed",
        #       n=np.arange(0, bootstrap_reps), archaic_genome=neanderthal_genomes),
        # expand(results_path + ldsc_path + '{archaic_genome}.heritability_estimates.txt', archaic_genome=neanderthal_genomes),
        # ## simulations and variant dating
        # # need to create 'data_path + AA_sample_ids.txt' that contains as many line as used in the actual analysis
        expand(simulations_path + "ALL_chromosomes_replicate_{n}_ancestry_proportions-{M}.{ext}",
               n=n_replicates, M=M, ext=['3.Q', 'fam']),
        expand(simulations_path + "AMR_introgression_deserts_replicate_{n}.bed", n=n_replicates),
        expand(simulations_path + "AMR_introgression_frequencies_and_rank_callable_windows_replicate_{n}.bed",
               n=n_replicates),
        expand(simulations_path + "AMR_novel_introgression_deserts_replicate_{n}.bed", n=n_replicates),
        expand(simulations_path + "AMR_novel_introgression_deserts_pvalues_replicate_{n}.bed", n=n_replicates),
        expand(simulations_path + "AMR_putatively_selected_neanderthal_segments" +
               "_replicate_{n}_window{windowsize}_s_{stepsize}.bed",
               n=n_replicates, windowsize=windowsize, stepsize=stepsize),
        expand(simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked" +
               "_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_pvalues.bed",
               n=n_replicates, windowsize=windowsize, stepsize=stepsize),
        expand(simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked" +
               "_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}_expectations.bed",
               n=n_replicates, windowsize=windowsize, stepsize=stepsize),
        expand(simulations_path + "neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked" +
               "_coverage_per_individual_and_per_window{windowsize}_s_{stepsize}.bed",
               windowsize=windowsize, n=n_replicates, stepsize=stepsize),
        expand(simulations_path + 'neanderthal_introgressed_segments_masked_denisovan_replicate_{n}_afr_masked.bed',
               n=n_replicates),
        expand(simulations_path +  "ALL_chromosomes_replicate_{n}_ancestry_proportions-{M}.{ext}",
               n=n_replicates, M=M, ext=['3.Q', 'fam']),
        # #expand(arg_dir + 'chr{chrom}_tmrca_summarized.tab', chrom=chromosomes)




