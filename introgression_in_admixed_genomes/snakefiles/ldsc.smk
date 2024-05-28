rule install_ldsc:
    output:
        'ldsc/ldsc.py'
    conda:
        '../envs/ldsc.yaml'
    shell:
        "git clone https://github.com/bulik/ldsc.git"

rule create_annotation_file_neg_selected:
    input:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments.bed"
    output:
        temp(results_path + ldsc_path + '{archaic_genome}_neg_selected.txt')
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes)
    shell:
        "tail -n+2 {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($4 < $11) print $1,$2 - 250000,$3 + 250000, \"neg_selected\"}}' "
        "> {output}"
        # pad by 100kb in both directions

rule create_annotation_file_pos_selected:
    input:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments.bed"
    output:
        temp(results_path + ldsc_path + '{archaic_genome}_pos_selected.txt')
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes)
    shell:
        "tail -n+2 {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($4 > $11) print $1,$2 - 250000,$3 + 250000, \"pos_selected\"}}' "
        "> {output}"
        # pad by 100kb in both directions

rule create_annotation_file_deserts:
    input:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts.bed"
    output:
        temp(results_path + ldsc_path + '{archaic_genome}_deserts.txt')
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes)
    shell:
        "tail -n+2 {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2 ,$3, \"desert\"}}' > {output}"

# def get_bootstrap(wildcards):
#     return wildcards.n

# rule create_annotation_file_not_selected:
#     input:
#         results_path + "ibdmix_{archaic_genome}/AMR_putatively_not_selected_control_neanderthal_segments_{n}.bed"
#     output:
#         temp(results_path + ldsc_path + '{archaic_genome}_not_selected_{n}.txt')
#     params:
#         bootstrap = get_bootstrap
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         n="[0-9]+"
#     shell:
#         "tail -n+2 {input} | awk -v rep=\"{params.bootstrap}\" -F "
#         "'\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1,$2 - 250000, $3 + 250000, \"not_selected_\"rep}}' > {output}"
        # pad by 100kb in both directions

# rule create_annotation_file_non_introgressed:
#     input:
#         results_path + "ibdmix_{archaic_genome}/AMR_putatively_non_introgressed_control_neanderthal_segments_{n}.bed"
#     output:
#         temp(results_path + ldsc_path + '{archaic_genome}_non_introgressed_{n}.txt')
#     params:
#         bootstrap = get_bootstrap
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         n="[0-9]+"
#     shell:
#         "tail -n+2 {input} | awk -v rep=\"{params.bootstrap}\" -F '\\t' "
#         "'BEGIN{{OFS=\"\\t\"}}{{print $1,$2 - 250000,$3 + 250000, \"non_introgressed_\"rep}}' > {output}"
        # pad by 100kb in both directions

#
# rule create_annotation_file_controls_deserts:
#     input:
#         results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_control_segments_{n}.bed"
#     output:
#         temp(results_path + ldsc_path + '{archaic_genome}_desert_control_{n}.txt')
#     params:
#         bootstrap = get_bootstrap
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         n="[0-9]+"
#     shell:
#         "tail -n+2 {input} | awk -v rep=\"{params.bootstrap}\" -F '\\t' "
#         "'BEGIN{{OFS=\"\\t\"}}{{print $1,$2 ,$3,\"desert_control_\"rep}}' > {output}"

rule download_hg38_to_hg19_chain:
    output:
        data_path + reference_path + "hg38ToHg19.over.chain.gz"
    params:
        data_path = data_path + reference_path,
        url = hg38_to_hg19_chain_url
    shell:
        "wget -q -P {params.data_path} {params.url}"

rule lift_annotation_file_to_hg19:
    input:
        bed=results_path + ldsc_path + '{archaic_genome}_{region_type}.txt',
        chain=data_path + reference_path + "hg38ToHg19.over.chain.gz"
    output:
        results_path + ldsc_path + '{archaic_genome}_{region_type}.bed'
    conda:
        "../envs/crossmap.yaml"
    wildcard_constraints:
        region_type='pos_selected|neg_selected|deserts'
    shell:
        'CrossMap.py bed --chromid s {input.chain} {input.bed} {output}'

# use rule lift_annotation_file_to_hg19 as lift_annotation_file_to_hg19_controls with:
#     input:
#         bed=results_path + ldsc_path + '{archaic_genome}_{region_type}_{n}.txt',
#         chain=data_path + reference_path + "hg38ToHg19.over.chain.gz"
#     output:
#         results_path + ldsc_path + '{archaic_genome}_{region_type}_{n}.bed'
#     wildcard_constraints:
#         region_type='not_selected|desert_control|non_introgressed',
#         n="[0-9]+"

rule get_eur_plink_file:
    output:
        expand(results_path + ldsc_path + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}{ext}",
               ext=[".bed", ".bim", ".fam"], chr=chromosomes)
    params:
        output_dir = results_path + ldsc_path,
        ldsc_plink_url = ldsc_plink_url,
        zipped = results_path + ldsc_path + ldsc_plink_url.split('/')[-1]
    shell:
        "wget -q -P {params.output_dir} {params.ldsc_plink_url}; tar xf {params.zipped} -C {params.output_dir}"


rule get_allele_frequencies_ldsc:
    output:
        expand(results_path + ldsc_path + ldsc_freq_url.split('/')[-1].replace('.tgz', '/') + "1000G.EUR.QC.{chr}.frq",
               chr=chromosomes)
    params:
        output_dir = results_path + ldsc_path,
        ldsc_freq_url = ldsc_freq_url,
        zipped=results_path + ldsc_path + ldsc_freq_url.split('/')[-1]
    shell:
        "wget -q -P {params.output_dir} {params.ldsc_freq_url}; tar xf {params.zipped} -C {params.output_dir}"

rule get_allele_weights_ldsc:
    output:
        expand(results_path + ldsc_path + ldsc_weights_url.split('/')[-1].replace('.tgz', '/') +
               "weights.hm3_noMHC.{chr}.l2.ldscore.gz", chr=chromosomes)
    params:
        output_dir = results_path + ldsc_path,
        ldsc_weights_url = ldsc_weights_url,
        zipped=results_path + ldsc_path + ldsc_weights_url.split('/')[-1]
    shell:
        "wget -q -P {params.output_dir} {params.ldsc_weights_url}; tar xf {params.zipped} -C {params.output_dir}"

rule get_baseline_annotation:
    output:
        annot = expand(results_path + ldsc_path + "baselineLD.{chr}.annot.gz", chr=chromosomes),
        aux = temp(expand(results_path + ldsc_path + "baselineLD.{chr}.{ext}", chr=chromosomes, ext=["l2.ldscore.gz",
                                                                                                     "l2.M", "l2.M_5_50",
                                                                                                     'log'])),
        tar = temp(results_path + ldsc_path + ldsc_baseline_model_url.split('/')[-1])
    params:
        data_dir = results_path + ldsc_path,
        url = ldsc_baseline_model_url,
        zipped = results_path + ldsc_path + ldsc_baseline_model_url.split('/')[-1]
    shell:
        "wget -q -P {params.data_dir} {params.url}; tar xf {params.zipped} -C {params.data_dir}"


rule make_ldsc_anno_file:
    input:
        pos_selected=results_path + ldsc_path + '{archaic_genome}_pos_selected.bed',
        neg_selected=results_path + ldsc_path + '{archaic_genome}_neg_selected.bed',
        deserts=results_path + ldsc_path + '{archaic_genome}_deserts.bed',
        baseline=results_path + ldsc_path + "baselineLD.{chr}.annot.gz"
    output:
        results_path + ldsc_path + "{archaic_genome}.{chr}.annot.gz"
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes),
        # region_type='selected|deserts'
    shell:
        "scripts/make_annot_ldsc.py --pos_selected {input.pos_selected} --neg_selected {input.neg_selected} "
        "--deserts {input.deserts} --baseline {input.baseline} --annot-file {output}"

# use rule make_ldsc_anno_file as make_ldsc_anno_file_controls with:
#     input:
#         coords=results_path + ldsc_path + '{archaic_genome}_{region_type}_{n}.bed',
#         baseline=results_path + ldsc_path + "baselineLD.{chr}.annot.gz"
#     output:
#         results_path + ldsc_path + "{archaic_genome}.{region_type}_regions_{n}.{chr}.annot.gz"
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         region_type='not_selected|desert_control|non_introgressed',
#         n="[0-9]+"

rule compute_ld_scores:
    input:
        plink=multiext(results_path + ldsc_path + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}", ".bed", ".bim", ".fam"),
        annot=results_path + ldsc_path + "{archaic_genome}.{chr}.annot.gz",
        LDSC='ldsc/ldsc.py'
    output:
        multiext(results_path + ldsc_path + "{archaic_genome}.{chr}.l2", ".M_5_50", ".ldscore.gz",
                 ".M")
    conda:
        "../envs/ldsc.yaml"
    params:
        plink_prefix=results_path + ldsc_path + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}",
        output_prefix=results_path + ldsc_path + "{archaic_genome}.{chr}"
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes),
        # region_type='selected|desert'
    resources:
        load = 5
    shell:
        "python ldsc/ldsc.py --l2 --bfile {params.plink_prefix} --ld-wind-cm 1 --annot {input.annot} --overlap-annot "
        "--out {params.output_prefix}"

# use rule compute_ld_scores as compute_ld_scores_controls with:
#     input:
#         plink=multiext(results_path + ldsc_path + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}", ".bed", ".bim", ".fam"),
#         annot=results_path + ldsc_path + "{archaic_genome}.{region_type}_regions_{n}.{chr}.annot.gz",
#         LDSC='ldsc/ldsc.py'
#     output:
#         multiext(results_path + ldsc_path + "{archaic_genome}.{region_type}_regions_{n}.{chr}.l2", ".M_5_50", ".ldscore.gz",
#                  ".M")
#     params:
#         plink_prefix=results_path + ldsc_path + "1000G_EUR_Phase3_plink/1000G.EUR.QC.{chr}",
#         output_prefix=results_path + ldsc_path + "{archaic_genome}.{region_type}_regions_{n}.{chr}"
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         region_type='not_selected|desert_control|non_introgressed',
#         n="[0-9]+"

checkpoint download_gwas_summary_stats:
    output:
        directory(results_path + ldsc_path + gwas_summary_stats_url.split('/')[-1].replace('.tgz', ''))
    params:
        output_dir = results_path + ldsc_path,
        zipped = results_path + ldsc_path + gwas_summary_stats_url.split('/')[-1],
        url = gwas_summary_stats_url
    shell:
        "wget -q -P {params.output_dir} {params.url}; tar xvf {params.zipped} -C {params.output_dir}"


def get_annotation_ldsc(wildcards):
    return [results_path + ldsc_path + f"{wildcards.archaic_genome}.{chrom}.annot.gz"
            for chrom in chromosomes]

def get_ldscores(wildcards):
    return [results_path + ldsc_path + f"{wildcards.archaic_genome}.{chrom}.l2.ldscore.gz"
            for chrom in chromosomes]


rule estimate_partitioned_heritability:
    input:
        weights=expand(results_path + ldsc_path + ldsc_weights_url.split('/')[-1].replace('.tgz', '/') +
                       "weights.hm3_noMHC.{chr}.l2.ldscore.gz", chr=chromosomes),
        sumstats=results_path + ldsc_path +  gwas_summary_stats_url.split('/')[-1].replace('.tgz', '/') + "{trait}.sumstats.gz",
        annot=get_annotation_ldsc,
        frequencies=expand(results_path + ldsc_path + ldsc_freq_url.split('/')[-1].replace('.tgz', '/')  +
                           "1000G.EUR.QC.{chr}.frq", chr=chromosomes),
        ldscores=get_ldscores
    output:
        multiext(results_path + ldsc_path +  "{archaic_genome}.{trait}", '.results', '.log')
    conda:
        "../envs/ldsc.yaml"
    params:
        weights_prefix=results_path + ldsc_path + ldsc_weights_url.split('/')[-1].replace('.tgz', '/') + "weights.hm3_noMHC.",
        freq_prefix=results_path + ldsc_path + ldsc_freq_url.split('/')[-1].replace('.tgz', '/') + "1000G.EUR.QC.",
        annot_prefix=results_path + ldsc_path + "{archaic_genome}.",
        output_prefix=results_path + ldsc_path + '{archaic_genome}.{trait}'
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes),
        # region_type='selected|desert'
    resources:
        load=5,
        mem='64GB'
    threads: 2
    shell:
        "python ldsc/ldsc.py --h2 {input.sumstats} --ref-ld-chr {params.annot_prefix} --w-ld-chr {params.weights_prefix} "
        "--overlap-annot --frqfile-chr {params.freq_prefix} --out {params.output_prefix}"

# def get_annotation_ldsc_controls(wildcards):
#     return [results_path + ldsc_path + f"{wildcards.archaic_genome}.{wildcards.region_type}_regions_{wildcards.n}.{chrom}.annot.gz"
#             for chrom in chromosomes]
#
#
# def get_ldscores_controls(wildcards):
#     return [results_path + ldsc_path + f"{wildcards.archaic_genome}.{wildcards.region_type}_regions_{wildcards.n}.{chrom}.l2.ldscore.gz"
#             for chrom in chromosomes]
#
# use rule estimate_partitioned_heritability as estimate_partitioned_heritability_controls with:
#     input:
#         weights=expand(results_path + ldsc_path + ldsc_weights_url.split('/')[-1].replace('.tgz', '/') +
#                        "weights.hm3_noMHC.{chr}.l2.ldscore.gz", chr=chromosomes),
#         sumstats=results_path + ldsc_path +  gwas_summary_stats_url.split('/')[-1].replace('.tgz', '/') + "{trait}.sumstats.gz",
#         annot=get_annotation_ldsc_controls,
#         frequencies=expand(results_path + ldsc_path + ldsc_freq_url.split('/')[-1].replace('.tgz', '/')  +
#                            "1000G.EUR.QC.{chr}.frq", chr=chromosomes),
#         ldscores=get_ldscores_controls
#     output:
#         multiext(results_path + ldsc_path +  "{archaic_genome}.{region_type}_regions_{n}_{trait}", '.results', '.log')
#     params:
#         weights_prefix=results_path + ldsc_path + ldsc_weights_url.split('/')[-1].replace('.tgz', '/') + "weights.hm3_noMHC.",
#         freq_prefix=results_path + ldsc_path + ldsc_freq_url.split('/')[-1].replace('.tgz', '/') + "1000G.EUR.QC.",
#         annot_prefix=results_path + ldsc_path + "{archaic_genome}.{region_type}_regions_{n}.",
#         output_prefix=results_path + ldsc_path + '{archaic_genome}.{region_type}_regions_{n}_{trait}'
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         region_type='not_selected|desert_control|non_introgressed',
#         n="[0-9]+"
#     resources:
#         load=5,
#         mem='64GB'
#     threads: 2


def aggregate_partitioned_heritability_estimates_input(wildcards):
    checkpoints.download_gwas_summary_stats.get(**wildcards)
    traits = glob_wildcards(results_path + ldsc_path +  gwas_summary_stats_url.split('/')[-1].replace('.tgz', '/') +
                            "{trait}.sumstats.gz").trait
    return expand(results_path + ldsc_path + "{archaic_genome}.{trait}.results",
                  trait=sorted(list(set(traits))), archaic_genome=wildcards.archaic_genome)

rule aggregate_partitioned_heritability_estimates:
    input:
        aggregate_partitioned_heritability_estimates_input
    output:
        results_path + ldsc_path + '{archaic_genome}.heritability_estimates.txt'
    wildcard_constraints:
        archaic_genome = "|".join(neanderthal_genomes),
        # region_type='selected|desert'
    shell:
        "echo {input} > {output}"

# def aggregate_partitioned_heritability_estimates_input_controls(wildcards):
#     checkpoints.download_gwas_summary_stats.get(**wildcards)
#     traits = glob_wildcards(results_path + ldsc_path +  gwas_summary_stats_url.split('/')[-1].replace('.tgz', '/') +
#                             "{trait}.sumstats.gz").trait
#     return expand(results_path + ldsc_path + "{archaic_genome}.{region_type}_regions_{n}_{trait}.results",
#                   trait=sorted(list(set(traits))), archaic_genome=wildcards.archaic_genome,
#                   region_type=wildcards.region_type, n=wildcards.n)
#
# rule aggregate_partitioned_heritability_estimates_controls:
#     input:
#         aggregate_partitioned_heritability_estimates_input_controls
#     output:
#         results_path + ldsc_path + '{archaic_genome}.{region_type}_{n}_regions_heritability_estimates.txt'
#     wildcard_constraints:
#         archaic_genome = "|".join(neanderthal_genomes),
#         region_type='not_selected|desert_control',#|non_introgressed',
#         n="[0-9]+"
#     shell:
#         "echo {input} > {output}"
