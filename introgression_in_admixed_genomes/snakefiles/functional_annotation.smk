import random
import string
import gzip
import pandas as pd

rule extract_genotypes_archaic_genomes:
    input:
        data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz"
    output:
        bed=data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.bed.gz",
        tbi=data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.bed.gz.tbi"
    params:
        bcftools = bcftools_path
    conda:
        "../envs/tabix.yaml"
    shell:
        "{params.bcftools} query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t[%GT]\\n' {input} | "
        "awk -F '\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1, $2 - 1, $2, $3, $4, $5}}' | "
        "sed -e 's/0\\/0/0/' -e 's/0\\/1/1/' -e 's/1\\/0/1/' -e 's/1\\/1/2/' | bgzip > {output.bed}; tabix {output.bed}"

def get_archaic_genotypes(wildcards):
    genomes = [data_path + "archaic_genomes/" + g + f'/vcf/chr{wildcards.chr}_mq25_mapab100.bed.gz' for g in archaic_genomes]
    return genomes

def get_archaic_genotype_indices(wildcards):
    indices = [data_path + "archaic_genomes/" + g + f'/vcf/chr{wildcards.chr}_mq25_mapab100.bed.gz.tbi' for g in archaic_genomes]
    return indices

rule create_vcfanno_config_file_modern:
    input:
        ancestral=data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_hg38_{chr}.bed.gz',
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        cds=data_path + reference_path + 'knownCDS.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','.bed.gz'),
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz', '.bed.gz'),
        phylop = data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed.gz',
        gnomad=data_path + reference_path + gnomad_url.split('/')[-1],
        clinvar=data_path + reference_path + clinvar_url.split('/')[-1],
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz', '.bed.gz'),
        archaic_genotypes=get_archaic_genotypes
    output:
        'config/vcfanno_chr{chr}.conf'
    params:
        anno_fields = anno_fields
    run:
        with open(output[0], 'w') as conf:
            for key, vals in anno_fields.items():
                try:
                   x = input[key]
                except:
                   continue
                conf.write("[[annotation]]\n")
                conf.write(f'file="{input[key]}"\n')
                for i, val in enumerate(vals):
                    field = val[0]
                    name= val[1]
                    op = val[2]
                    if input[key].endswith('.bed.gz') and i == 0:
                        fields = f"columns=[{field}"
                    elif input[key].endswith('.vcf.gz') and i == 0:
                        fields = f'fields=["{field}"'
                    elif input[key].endswith('.bed.gz') and i > 0:
                        fields += f", {field}"
                    elif input[key].endswith('.vcf.gz') and i > 0:
                        fields += f', "{field}"'
                    if i == 0:
                        names = f'names=["{name}"'
                        ops = f'ops=["{op}"'
                    elif i > 0:
                        names += f', "{name}"'
                        ops += f', "{op}"'
                fields += "]"
                names += "]"
                ops += "]"
                conf.write(f"{fields}\n")
                conf.write(f"{names}\n")
                conf.write(f"{ops}\n")
                conf.write("\n")
            for archaic_genome in input.archaic_genotypes:
                conf.write("[[annotation]]\n")
                conf.write(f'file="{archaic_genome}"\n')
                conf.write(f'columns=[4, 5, 6]\n')
                conf.write(f'names=["{archaic_genome.split("/")[-3]}_REF", "{archaic_genome.split("/")[-3]}_ALT", '
                                    f'"{archaic_genome.split("/")[-3]}_GT"]\n')
                conf.write('ops=["self", "self", "self"]\n')
                conf.write('\n')
        conf.close()


rule create_vcfanno_config_file_modern_1kgp_phased:
    input:
        ancestral=data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_hg38_{chr}.bed.gz',
    output:
        temp('config/vcfanno_1kgp_chr{chr}.conf')
    params:
        anno_fields = anno_fields
    run:
        with open(output[0], 'w') as conf:
            key = 'ancestral'
            vals = anno_fields[key]
            conf.write("[[annotation]]\n")
            conf.write(f'file="{input[key]}"\n')
            for i, val in enumerate(vals):
                field = val[0]
                name= val[1]
                op = val[2]
                if i == 0:
                    names = f'names=["{name}"'
                    ops = f'ops=["{op}"'
                    fields = f"columns=[{field}"
                elif i > 0:
                    fields += f", {field}"
                    names += f', "{name}"'
                    ops += f', "{op}"'
            fields += "]"
            names += "]"
            ops += "]"
            conf.write(f"{fields}\n")
            conf.write(f"{names}\n")
            conf.write(f"{ops}\n")
            conf.write("\n")
        conf.close()

rule annotate_modern_vcf:
    input:
        target_vcf=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        annotation_config='config/vcfanno_chr{chr}.conf',
        genes=multiext(data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz', ''),
                       '.bed.gz', '.bed.gz.tbi'),
        cds=multiext(data_path + reference_path + 'knownCDS', '.bed.gz', '.bed.gz.tbi'),
        gerp=multiext(data_path + reference_path + gerp_url.split('/')[-1].replace('.bw', ''), '.bed.gz', '.bed.gz.tbi'),
        phastcons=multiext(data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz', ''),
                           '.bed.gz', '.bed.gz.tbi'),
        phylop=multiext(data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz', ''),
                           '.bed.gz', '.bed.gz.tbi'),
        recomb=multiext(data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', ''),
                           '.bed.gz', '.bed.gz.tbi'),
        bstat=multiext(data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38'), '.bed.gz', '.bed.gz.tbi'),
        gnomad=data_path + reference_path + gnomad_url.split('/')[-1],
        gnomad_tbi=data_path + reference_path + gnomad_url.split('/')[-1] + '.tbi',
        clinvar=data_path + reference_path + clinvar_url.split('/')[-1],
        clinvar_tbi=data_path + reference_path + clinvar_url.split('/')[-1] + '.tbi',
        encode=multiext(data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz', ''),
                        '.bed.gz', '.bed.gz.tbi'),
        archaic_genotypes=get_archaic_genotypes,
        ancestral=data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_hg38_{chr}.bed.gz',
        archaic_genotype_indices=get_archaic_genotype_indices
    output:
        temp(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated.vcf")
    conda:
        "../envs/vcfanno.yaml"
    threads: 8
    shell:
        "vcfanno -p {threads} {input.annotation_config} {input.target_vcf} > {output}"


rule compress_and_index_annotated_modern_vcf:
    input:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated.vcf"
    output:
        vcf=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated.vcf.gz",
        tbi=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated.vcf.gz.tbi"
    conda:
        '../envs/tabix.yaml'
    shell:
        "bgzip {input}; tabix {output.vcf}"

rule drop_genotypes:
    input:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated.vcf.gz"
    output:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz"
    params:
        bcftools = bcftools_path
    shell:
        "{params.bcftools} view -G -Oz -o {output} {input}"

rule index_vcf_without_genotypes:
    input:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz"
    output:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz.tbi"
    conda:
        '../envs/tabix.yaml'
    shell:
        "tabix {input}"

rule generate_config_file_interval_annotation:
    input:
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','') + '.bed.gz',
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        phylop=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed.gz',
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        cds=data_path + reference_path + 'knownCDS.bed.gz',
        boosting_complete_eur=data_path + reference_path + "EUR_complete_boosting_scores.bed.gz",
        boosting_complete_eur_recent=data_path + reference_path + "EUR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eur_ancient=data_path + reference_path + "EUR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eur=data_path + reference_path + "EUR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eur_recent=data_path + reference_path + "EUR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eur_ancient=data_path + reference_path + "EUR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_afr=data_path + reference_path + "AFR_complete_boosting_scores.bed.gz",
        boosting_complete_afr_recent=data_path + reference_path + "AFR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_afr_ancient=data_path + reference_path + "AFR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_afr=data_path + reference_path + "AFR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_afr_recent=data_path + reference_path + "AFR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_afr_ancient=data_path + reference_path + "AFR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_eas=data_path + reference_path + "EAS_complete_boosting_scores.bed.gz",
        boosting_complete_eas_recent=data_path + reference_path + "EAS_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eas_ancient=data_path + reference_path + "EAS_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eas=data_path + reference_path + "EAS_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eas_recent=data_path + reference_path + "EAS_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eas_ancient=data_path + reference_path + "EAS_incomplete_ancient_boosting_scores.bed.gz",
        pi=results_path + "nucleotide_diversity_per_10kb.bed.gz"
    output:
        temp("config/config_interval_anno.txt")
    params:
        anno_fields = anno_fields
    run:
        with open(output[0], 'w') as conf:
            for key, vals in anno_fields.items():
                try:
                    conf.write(f'{input[key]}\t')
                except AttributeError:
                    conf.write(f'{input[key]}\t')
                bed = gzip.open(input[key], 'r')
                n_fields = len(bed.readline().strip().split(b'\t'))
                bed.close()
                fields = ['chrom_b', 'start_b', 'end_b']
                targets = []
                field = int(vals[0][0])
                name = vals[0][1]
                data_type = vals[0][3]
                if field == len(fields) + 1:
                    fields.append(name)
                    targets.append(name)
                else:
                    while field > len(fields) + 1:
                        fields.append("".join(random.choices(string.ascii_lowercase, k=8)))
                    fields.append(name)
                    targets.append(name)
                while len(fields) < n_fields:
                    fields.append("".join(random.choices(string.ascii_lowercase,k=8)))
                conf.write(','.join(fields) + '\t')
                conf.write(','.join(targets) + '\t')
                conf.write(data_type)
                conf.write('\n')
        conf.close()

rule merge_known_deserts:
    input:
        data_path + "known_introgression_deserts_hg38.bed"
    output:
        results_path + "known_introgression_deserts_hg38_merged.bed"
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} merge -i {input} > {output}"


rule calculate_idat:
    input:
        flare_0=results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase0.{archaic_genome}.bed",
        flare_1=results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase1.{archaic_genome}.bed",
        genomefile=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        results_path + 'ibdmix_{archaic_genome}/iDAT_scores_{chr}.bed',
    params:
        chrom="{chr}",
        total=n_sites_idat,
        dist_step_size=dist_step_size,
        min_cov_dat=min_cov_dat
    threads: 8
    shell:
        "scripts/calculate_iDAT_scores.py idat_scores -l {input.flare_0} {input.flare_1} -g {input.genomefile} "
        "-o {output} -c {params.chrom} -n {params.total} -d {params.dist_step_size} -m {params.min_cov_dat} "
        "--threads {threads}"

rule calculate_standardized_idat_windows:
    input:
        idat=[results_path + 'ibdmix_{archaic_genome}/' + f'iDAT_scores_{chrom}.bed' for chrom in chromosomes],
        genomefile=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        results_path + 'ibdmix_{archaic_genome}/standardized_iDAT_scores.bed'
    params:
        windowsize=windowsize_idat,
        stride=stepsize_idat
    shell:
        "scripts/calculate_iDAT_scores.py standardized_idat_scores -i {input.idat} -g {input.genomefile} -o {output} "
        "-w {params.windowsize} -s {params.stride}"

rule calculate_standardized_idat_selected_regions:
    input:
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}." +
             "{archaic_genome}.bed" for chrom in chromosomes for phase in [0, 1]],
        regions=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments.bed",
        idat=[results_path + 'ibdmix_{archaic_genome}/' + f'iDAT_scores_{chrom}.bed' for chrom in chromosomes]
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed"
    params:
        lai_pattern=results_path + flare_output + "african_american_and_ref_individuals_chr1.anc_per_pos.phase0.{archaic_genome}.bed",
        dist_step_size=dist_step_size,
        min_cov_dat=min_cov_dat
    threads: 16
    shell:
        "scripts/calculate_iDAT_scores.py regions -l {params.lai_pattern} -r {input.regions} --idat {input.idat} "
        "-o {output} -d {params.dist_step_size} -m {params.min_cov_dat} --threads {threads}"

use rule calculate_standardized_idat_selected_regions as calculate_standardized_idat_regions with:
    input:
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}." +
             "{archaic_genome}.bed" for chrom in chromosomes for phase in [0, 1]],
        regions=results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}.bed",
        idat=[results_path + 'ibdmix_{archaic_genome}/' + f'iDAT_scores_{chrom}.bed' for chrom in chromosomes]
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed"
    params:
        lai_pattern=results_path + flare_output + "african_american_and_ref_individuals_chr1.anc_per_pos.phase0.{archaic_genome}.bed",
        dist_step_size=dist_step_size,
        min_cov_dat=min_cov_dat
    wildcard_constraints:
        segment_type='not_selected_control|non_introgressed_control',
        n="[0-9][0-9]?"
    threads: 16

use rule calculate_standardized_idat_selected_regions as calculate_standardized_idat_deserts with:
    input:
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}." +
             "{archaic_genome}.bed" for chrom in chromosomes for phase in [0, 1]],
        regions=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts.bed",
        idat=[results_path + 'ibdmix_{archaic_genome}/' + f'iDAT_scores_{chrom}.bed' for chrom in chromosomes]
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed"
    params:
        lai_pattern=results_path + flare_output + "african_american_and_ref_individuals_chr1.anc_per_pos.phase0.{archaic_genome}.bed",
        dist_step_size=dist_step_size,
        min_cov_dat=min_cov_dat
    threads: 16

use rule calculate_standardized_idat_selected_regions as calculate_standardized_idat_deserts_control_regions with:
    input:
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}." +
             "{archaic_genome}.bed" for chrom in chromosomes for phase in [0, 1]],
        regions=results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}.bed",
        idat=[results_path + 'ibdmix_{archaic_genome}/' +  f'iDAT_scores_{chrom}.bed' for chrom in chromosomes]
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed"
    params:
        lai_pattern=results_path + flare_output + "african_american_and_ref_individuals_chr1.anc_per_pos.phase0.{archaic_genome}.bed",
        dist_step_size=dist_step_size,
        min_cov_dat=min_cov_dat
    threads: 16
    wildcard_constraints:
        n = "[0-9][0-9]?"

use rule calculate_standardized_idat_selected_regions as calculate_standardized_idat_known_deserts with:
    input:
        lai=[results_path + flare_output + f"african_american_and_ref_individuals_chr{chrom}.anc_per_pos.phase{phase}." +
             "{archaic_genome}.bed" for chrom in chromosomes for phase in [0, 1]],
        regions=results_path + "known_introgression_deserts_hg38_merged.bed",
        idat=[results_path + 'ibdmix_{archaic_genome}/' + f'iDAT_scores_{chrom}.bed' for chrom in chromosomes]
    output:
        results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_iDAT_annotated.bed"
    params:
        lai_pattern=results_path + flare_output + "african_american_and_ref_individuals_chr1.anc_per_pos.phase0.{archaic_genome}.bed",
        dist_step_size=dist_step_size,
        min_cov_dat=min_cov_dat
    threads: 16


rule calculate_pi:
    input:
        data_path + "target_individuals_chr{chr}.vcf.gz"
    output:
        data_path + "target_individuals_chr{chr}.windowed.pi"
    conda:
        "../envs/vcftools.yaml"
    params:
        out=data_path + "target_individuals_chr{chr}."
    shell:
        "vcftools --gzvcf {input} --window-pi 10000 --out {params.out}"

rule aggregate_pi:
    input:
        expand(data_path + "target_individuals_chr{chr}.windowed.pi", chr=chromosomes)
    output:
        results_path + "nucleotide_diversity_per_10kb.bed.gz"
    conda:
        '../envs/tabix.yaml'
    shell:
        "tail -n+2 {input} | sort -k1,1 -k2,2n | bgzip > {output}"


# rule calculate_pi_selected:
#     input:
#         vcf=data_path + 'target_individuals_chr{chr}.vcf.gz',
#         regions=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
#         outliers=results_path + 'ibdmix_{archaic_genome}/admixed_outlier_iids.txt'
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_chr{chr}.sites.pi"
#     conda:
#         "../envs/vcftools.yaml"
#     params:
#         chrom=f'chr{chr}',
#         output=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_chr{chr}"
#     shell:
#         "vcftools --gzvcf {input.vcf} --site-pi --remove {input.outliers} "
#         "--bed <(grep -w {params.chrom} {input.regions} | cut-f1-3) --output {params.output}"
#
# use rule calculate_pi_selected as calculate_pi_not_selected with:
#     input:
#         vcf=data_path + 'target_individuals_chr{chr}.vcf.gz',
#         regions=results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
#         outliers=results_path + 'ibdmix_{archaic_genome}/admixed_outlier_iids.txt'
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_chr{chr}_{n}.sites.pi"
#     conda:
#         "../envs/vcftools.yaml"
#     wildcard_constraints:
#         chr="|".join([str(chrom) for chrom in chromosomes]),
#         n="[0-9]+"
#     params:
#         chrom=f'chr{chr}',
#         output=results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_chr{chr}_{n}"
#
# use rule calculate_pi_selected as calculate_pi_deserts with:
#     input:
#         vcf=data_path + 'target_individuals_chr{chr}.vcf.gz',
#         regions=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
#         outliers=results_path + 'ibdmix_{archaic_genome}/admixed_outlier_iids.txt'
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_chr{chr}.sites.pi"
#     conda:
#         "../envs/vcftools.yaml"
#     params:
#         chrom=f'chr{chr}',
#         output=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_desert_chr{chr}"
#
# use rule calculate_pi_selected as calculate_pi_deserts_control with:
#     input:
#         vcf=data_path + 'target_individuals_chr{chr}.vcf.gz',
#         regions=results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
#         outliers=results_path + 'ibdmix_{archaic_genome}/admixed_outlier_iids.txt'
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_chr{chr}_{n}.sites.pi"
#     conda:
#         "../envs/vcftools.yaml"
#     wildcard_constraints:
#         chr = "|".join([str(chrom) for chrom in chromosomes]),
#         n="[0-9]+"
#     params:
#         chrom=f'chr{chr}',
#         output=results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_chr{chr}_{n}"
#
# def get_pi_stats_files_selected(wildcards):
#     return [results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_putatively_selected_neanderthal_segments_chr{chrom}.sites.pi"
#             for chrom in chromosomes]
#
# def get_pi_stats_files_selected_control(wildcards):
#     return [results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_putatively_{wildcards.segment_type}_neanderthal_segments_chr{chrom}_{wildcards.n}.sites.pi"
#             for chrom in chromosomes]
#
# def get_pi_stats_files_deserts(wildcards):
#     return [results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_novel_introgression_deserts_chr{chrom}.sites.pi"
#             for chrom in chromosomes]
#
# def get_pi_stats_files_deserts_control(wildcards):
#     return [results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_introgression_deserts_new_control_segments_chr{chrom}_{wildcards.n}.sites.pi"
#             for chrom in chromosomes]
#
# rule aggregate_pi_statistic_selected:
#     input:
#         get_pi_stats_files_selected
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_sites_pi.bed.gz"
#     conda:
#         '../envs/tabix.yaml'
#     shell:
#         "tail -n+2 {input} | sort -k1,1 -k2,2n | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1, $2-1, $2, $3}}' | bgzip > {output}"
#
# use rule aggregate_pi_statistic_selected as aggregate_pi_statistic_selected_control with:
#     input:
#         get_pi_stats_files_selected_control
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_sites_pi.bed.gz"
#     wildcard_constraints:
#         segment_type = 'not_selected_control|non_introgressed_control',
#         n="[0-9]+"
#     conda:
#         '../envs/tabix.yaml'
#
# use rule aggregate_pi_statistic_selected as aggregate_pi_statistic_deserts with:
#     input:
#         get_pi_stats_files_deserts
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_sites_pi.bed.gz"
#     conda:
#         '../envs/tabix.yaml'
#
# use rule aggregate_pi_statistic_selected as aggregate_pi_statistic_deserts_control with:
#     input:
#         get_pi_stats_files_deserts_control
#     output:
#         results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_sites_pi.bed.gz"
#     conda:
#         '../envs/tabix.yaml'
#     wildcard_constraints:
#         n="[0-9]+"

# only using bed files here, I cannot use annotated VCF of modern humans as they only contain a subset of the
# data points, i.e., only infos about sites that are variable, while using the raw files better allows characterizing
# the regions as a whole
rule annotate_selected_introgressed_intervals_basic_stats:
    input:
        conf = "config/config_interval_anno.txt",
        nea = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','') + '.bed.gz',
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        phylop=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed.gz',
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        boosting_complete_eur=data_path + reference_path + "EUR_complete_boosting_scores.bed.gz",
        boosting_complete_eur_recent=data_path + reference_path + "EUR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eur_ancient=data_path + reference_path + "EUR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eur=data_path + reference_path + "EUR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eur_recent=data_path + reference_path + "EUR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eur_ancient=data_path + reference_path + "EUR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_afr=data_path + reference_path + "AFR_complete_boosting_scores.bed.gz",
        boosting_complete_afr_recent=data_path + reference_path + "AFR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_afr_ancient=data_path + reference_path + "AFR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_afr=data_path + reference_path + "AFR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_afr_recent=data_path + reference_path + "AFR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_afr_ancient=data_path + reference_path + "AFR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_eas=data_path + reference_path + "EAS_complete_boosting_scores.bed.gz",
        boosting_complete_eas_recent=data_path + reference_path + "EAS_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eas_ancient=data_path + reference_path + "EAS_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eas=data_path + reference_path + "EAS_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eas_recent=data_path + reference_path + "EAS_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eas_ancient=data_path + reference_path + "EAS_incomplete_ancient_boosting_scores.bed.gz",
        pi=results_path + "nucleotide_diversity_per_10kb.bed.gz"
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_annotated_tmp.bed")
    resources:
        load = 10,
        mem_mb=64 * 1000
    shell:
        'scripts/annotate_segments_basic_stats.py -i {input.nea} -c {input.conf} --pi {input.pi} -o {output}'

use rule annotate_selected_introgressed_intervals_basic_stats as annotate_introgressed_intervals_basic_stats with:
    input:
        conf = "config/config_interval_anno.txt",
        nea = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','') + '.bed.gz',
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        phylop=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38')  + '.bed.gz',
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        boosting_complete_eur=data_path + reference_path + "EUR_complete_boosting_scores.bed.gz",
        boosting_complete_eur_recent=data_path + reference_path + "EUR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eur_ancient=data_path + reference_path + "EUR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eur=data_path + reference_path + "EUR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eur_recent=data_path + reference_path + "EUR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eur_ancient=data_path + reference_path + "EUR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_afr=data_path + reference_path + "AFR_complete_boosting_scores.bed.gz",
        boosting_complete_afr_recent=data_path + reference_path + "AFR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_afr_ancient=data_path + reference_path + "AFR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_afr=data_path + reference_path + "AFR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_afr_recent=data_path + reference_path + "AFR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_afr_ancient=data_path + reference_path + "AFR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_eas=data_path + reference_path + "EAS_complete_boosting_scores.bed.gz",
        boosting_complete_eas_recent=data_path + reference_path + "EAS_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eas_ancient=data_path + reference_path + "EAS_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eas=data_path + reference_path + "EAS_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eas_recent=data_path + reference_path + "EAS_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eas_ancient=data_path + reference_path + "EAS_incomplete_ancient_boosting_scores.bed.gz",
        pi=results_path + "nucleotide_diversity_per_10kb.bed.gz"
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_annotated_tmp.bed")
    resources:
        load = 10,
        mem_mb=64 * 1000
    wildcard_constraints:
        segment_type = 'not_selected_control|non_introgressed_control',

rule annotate_introgression_deserts_basic_stats:
    input:
        conf = "config/config_interval_anno.txt",
        nea = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','') + '.bed.gz',
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        phylop=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38')  + '.bed.gz',
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        boosting_complete_eur=data_path + reference_path + "EUR_complete_boosting_scores.bed.gz",
        boosting_complete_eur_recent=data_path + reference_path + "EUR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eur_ancient=data_path + reference_path + "EUR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eur=data_path + reference_path + "EUR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eur_recent=data_path + reference_path + "EUR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eur_ancient=data_path + reference_path + "EUR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_afr=data_path + reference_path + "AFR_complete_boosting_scores.bed.gz",
        boosting_complete_afr_recent=data_path + reference_path + "AFR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_afr_ancient=data_path + reference_path + "AFR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_afr=data_path + reference_path + "AFR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_afr_recent=data_path + reference_path + "AFR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_afr_ancient=data_path + reference_path + "AFR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_eas=data_path + reference_path + "EAS_complete_boosting_scores.bed.gz",
        boosting_complete_eas_recent=data_path + reference_path + "EAS_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eas_ancient=data_path + reference_path + "EAS_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eas=data_path + reference_path + "EAS_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eas_recent=data_path + reference_path + "EAS_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eas_ancient=data_path + reference_path + "EAS_incomplete_ancient_boosting_scores.bed.gz",
        pi=results_path + "nucleotide_diversity_per_10kb.bed.gz"
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_annotated_tmp.bed")
    resources:
        load = 10,
        mem_mb=64 * 1000
    shell:
        'scripts/annotate_segments_basic_stats.py -i {input.nea} -c {input.conf} -o {output} --quantitative'

use rule annotate_introgression_deserts_basic_stats as annotate_introgression_desert_control_segments_basic_stats with:
    input:
        conf = "config/config_interval_anno.txt",
        nea = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','') + '.bed.gz',
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        phylop=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38')  + '.bed.gz',
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        boosting_complete_eur=data_path + reference_path + "EUR_complete_boosting_scores.bed.gz",
        boosting_complete_eur_recent=data_path + reference_path + "EUR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eur_ancient=data_path + reference_path + "EUR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eur=data_path + reference_path + "EUR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eur_recent=data_path + reference_path + "EUR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eur_ancient=data_path + reference_path + "EUR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_afr=data_path + reference_path + "AFR_complete_boosting_scores.bed.gz",
        boosting_complete_afr_recent=data_path + reference_path + "AFR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_afr_ancient=data_path + reference_path + "AFR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_afr=data_path + reference_path + "AFR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_afr_recent=data_path + reference_path + "AFR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_afr_ancient=data_path + reference_path + "AFR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_eas=data_path + reference_path + "EAS_complete_boosting_scores.bed.gz",
        boosting_complete_eas_recent=data_path + reference_path + "EAS_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eas_ancient=data_path + reference_path + "EAS_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eas=data_path + reference_path + "EAS_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eas_recent=data_path + reference_path + "EAS_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eas_ancient=data_path + reference_path + "EAS_incomplete_ancient_boosting_scores.bed.gz",
        pi=results_path + "nucleotide_diversity_per_10kb.bed.gz"
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_annotated_tmp.bed")
    resources:
        load = 10,
        mem_mb=64 * 1000

use rule annotate_introgression_deserts_basic_stats as annotate_introgression_known_deserts_basic_stats with:
    input:
        conf = "config/config_interval_anno.txt",
        nea = results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_iDAT_annotated.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        gerp=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','') + '.bed.gz',
        phastcons=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        phylop=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        recomb=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        bstat=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38')  + '.bed.gz',
        encode=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','') + '.bed.gz',
        boosting_complete_eur=data_path + reference_path + "EUR_complete_boosting_scores.bed.gz",
        boosting_complete_eur_recent=data_path + reference_path + "EUR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eur_ancient=data_path + reference_path + "EUR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eur=data_path + reference_path + "EUR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eur_recent=data_path + reference_path + "EUR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eur_ancient=data_path + reference_path + "EUR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_afr=data_path + reference_path + "AFR_complete_boosting_scores.bed.gz",
        boosting_complete_afr_recent=data_path + reference_path + "AFR_complete_recent_boosting_scores.bed.gz",
        boosting_complete_afr_ancient=data_path + reference_path + "AFR_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_afr=data_path + reference_path + "AFR_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_afr_recent=data_path + reference_path + "AFR_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_afr_ancient=data_path + reference_path + "AFR_incomplete_ancient_boosting_scores.bed.gz",
        boosting_complete_eas=data_path + reference_path + "EAS_complete_boosting_scores.bed.gz",
        boosting_complete_eas_recent=data_path + reference_path + "EAS_complete_recent_boosting_scores.bed.gz",
        boosting_complete_eas_ancient=data_path + reference_path + "EAS_complete_ancient_boosting_scores.bed.gz",
        boosting_incomplete_eas=data_path + reference_path + "EAS_incomplete_boosting_scores.bed.gz",
        boosting_incomplete_eas_recent=data_path + reference_path + "EAS_incomplete_recent_boosting_scores.bed.gz",
        boosting_incomplete_eas_ancient=data_path + reference_path + "EAS_incomplete_ancient_boosting_scores.bed.gz",
        pi=results_path + "nucleotide_diversity_per_10kb.bed.gz"
    output:
        temp(results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_annotated_tmp.bed")
    resources:
        load = 10,
        mem_mb=64 * 1000

rule get_allele_frequencies_reference:
    input:
        plink=multiext(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.updated_ids", '.bed', '.bim', '.fam'),
        samples=data_path + "{superpopulation}_sample_ids.txt",
        fam=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.updated_ids.fam"
    output:
        data_path + "{superpopulation}_chr{chr}.afreq"
    conda:
        "../envs/plink.yaml"
    params:
        input_base = data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.updated_ids",
        output_base = data_path + "{superpopulation}_chr{chr}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    wildcard_constraints:
        superpopulation='AFR|EUR|EAS'
    shell:
        "plink2 --bfile {params.input_base} --keep <( grep -f {input.samples} {input.fam} | cut -f1,2) "
        "--freq --out {params.output_base} --threads {threads} --memory {resources.mem_mb}"


use rule get_allele_frequencies_reference as get_allele_frequencies_target with:
    input:
        plink=multiext(data_path + 'target_individuals_chr{chr}', '.bed', '.bim', '.fam'),
        samples=data_path + "AA_sample_ids.txt",
        fam=data_path + 'target_individuals_chr{chr}.fam'
    output:
        data_path + "AA_chr{chr}.afreq"
    conda:
        "../envs/plink.yaml"
    params:
        input_base = data_path + 'target_individuals_chr{chr}',
        output_base = data_path + "AA_chr{chr}"
    threads: 8
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000

def get_number_of_cols_selected(wildcards):
    columns = pd.read_csv(results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
                          header=0,nrows=1,sep='\t').columns.values.tolist()
    return len(columns) + 2


rule annotate_selected_introgressed_intervals_modern_variants:
    input:
        nea = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
        vcfs=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz", chr=chromosomes),
        tbi=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz.tbi", chr=chromosomes)
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_annotated_tmp2.bed")
    params:
        bedtools = bedtools_path,
        ncols=get_number_of_cols_selected
    shell:
        '{params.bedtools} intersect -wao -a <(tail -n+2 {input.nea} | sort -k1,1 -k2,2n) -b {input.vcfs} -sorted | cut -f1-3,{params.ncols}- > {output}'

def get_number_of_cols(wildcards):
    columns = pd.read_csv(results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_putatively_{wildcards.segment_type}_neanderthal_segments_iDAT_{wildcards.n}_annotated.bed",
                          header=0,nrows=1,sep='\t').columns.values.tolist()
    return len(columns) + 2

use rule annotate_selected_introgressed_intervals_modern_variants as annotate_introgressed_intervals_modern_variants with:
    input:
        nea = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
        vcfs=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz",chr=chromosomes),
        tbi=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz.tbi",chr=chromosomes)
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_annotated_tmp2.bed")
    params:
        bedtools = bedtools_path,
        ncols=get_number_of_cols
    wildcard_constraints:
        segment_type = 'not_selected_control|non_introgressed_control',

def get_number_of_cols_deserts(wildcards):
    columns = pd.read_csv(results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
                          header=0,nrows=1,sep='\t').columns.values.tolist()
    return len(columns) + 2

use rule annotate_selected_introgressed_intervals_modern_variants as annotate_introgression_deserts_modern_variants with:
    input:
        nea = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
        vcfs=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz",chr=chromosomes),
        tbi=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz.tbi",chr=chromosomes)
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_annotated_tmp2.bed")
    params:
        bedtools = bedtools_path,
        ncols=get_number_of_cols_deserts

def get_number_of_cols_desert_controls(wildcards):
    columns = pd.read_csv(results_path + f"ibdmix_{wildcards.archaic_genome}/AMR_introgression_deserts_new_control_segments_{wildcards.n}_iDAT_annotated.bed",
                          header=0,nrows=1,sep='\t').columns.values.tolist()
    return len(columns) + 2

use rule annotate_selected_introgressed_intervals_modern_variants as annotate_introgression_desert_control_segments_modern_variants with:
    input:
        nea = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        vcfs=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz",chr=chromosomes),
        tbi=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz.tbi",chr=chromosomes)
    output:
        temp(results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_annotated_tmp2.bed")
    params:
        bedtools = bedtools_path,
        ncols=get_number_of_cols_desert_controls

def get_number_of_cols_known_deserts(wildcards):
    columns = pd.read_csv(results_path + f"ibdmix_{wildcards.archaic_genome}/known_introgression_deserts_hg38_iDAT_annotated.bed",
                          header=0,nrows=1,sep='\t').columns.values.tolist()
    return len(columns) + 2

use rule annotate_selected_introgressed_intervals_modern_variants as annotate_introgression_known_deserts_modern_variants with:
    input:
        nea = results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_iDAT_annotated.bed",
        vcfs=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz",chr=chromosomes),
        tbi=expand(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz.tbi",chr=chromosomes)
    output:
        temp(results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_annotated_tmp2.bed")
    params:
        bedtools = bedtools_path,
        ncols=get_number_of_cols_known_deserts

rule get_vcf_cols_as_header:
    input:
        data_path + "1000G_phase3/REF.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.annotated_no_genotypes.vcf.gz"
    output:
        "header_annotated_vcf.txt"
    params:
        bcftools = bcftools_path
    shell:
        "{params.bcftools} view -h {input} | tail -n1 > {output}"

rule aggregate_modern_variant_annotations_per_selected_segment:
    input:
        modern_var=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_annotated_tmp2.bed",
        basic_stat=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_annotated_tmp.bed",
        vcf_cols='header_annotated_vcf.txt',
        aa_af=expand(data_path + "AA_chr{chr}.afreq", chr=chromosomes),
        afr_af=expand(data_path + "AFR_chr{chr}.afreq", chr=chromosomes),
        eur_af=expand(data_path + "EUR_chr{chr}.afreq", chr=chromosomes),
        eas_af=expand(data_path + "EAS_chr{chr}.afreq", chr=chromosomes)
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_annotated.bed"
    params:
        archaic_genomes = " ".join(archaic_genomes)
    resources:
            mem_mb=64 * 1000,
            load=10
    shell:
        'scripts/aggregate_annotation_of_modern_variants_putatively_selected_introgressed_regions.py '
        '-i {input.modern_var} -b {input.basic_stat} --header_modern_stats {input.vcf_cols} --aa_af {input.aa_af} '
        '--afr_af {input.afr_af} --eur_af {input.eur_af} --eas_af {input.eas_af} '
        '--archaic_genomes {params.archaic_genomes} --output {output}'

use rule aggregate_modern_variant_annotations_per_selected_segment as aggregate_modern_variant_annotations_per_segment with:
    input:
        modern_var=results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_annotated_tmp2.bed",
        basic_stat=results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_annotated_tmp.bed",
        vcf_cols='header_annotated_vcf.txt',
        aa_af=expand(data_path + "AA_chr{chr}.afreq", chr=chromosomes),
        afr_af=expand(data_path + "AFR_chr{chr}.afreq", chr=chromosomes),
        eur_af=expand(data_path + "EUR_chr{chr}.afreq", chr=chromosomes),
        eas_af=expand(data_path + "EAS_chr{chr}.afreq", chr=chromosomes)
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_annotated.bed"
    params:
        archaic_genomes = " ".join(archaic_genomes)
    resources:
        mem_mb=64 * 1000,
        load=10
    wildcard_constraints:
        segment_type = "not_selected_control|non_introgressed_control",
        n="[0-9]+"


use rule aggregate_modern_variant_annotations_per_selected_segment as aggregate_modern_variant_annotations_per_desert with:
    input:
        modern_var=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_annotated_tmp2.bed",
        basic_stat=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_annotated_tmp.bed",
        vcf_cols='header_annotated_vcf.txt',
        aa_af=expand(data_path + "AA_chr{chr}.afreq", chr=chromosomes),
        afr_af=expand(data_path + "AFR_chr{chr}.afreq", chr=chromosomes),
        eur_af=expand(data_path + "EUR_chr{chr}.afreq", chr=chromosomes),
        eas_af=expand(data_path + "EAS_chr{chr}.afreq", chr=chromosomes)
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_annotated.bed"
    params:
        archaic_genomes = " ".join(archaic_genomes)
    resources:
        mem_mb=64 * 1000,
        load=10

use rule aggregate_modern_variant_annotations_per_selected_segment as aggregate_modern_variant_annotations_per_desert_control_segment with:
    input:
        modern_var=results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_annotated_tmp2.bed",
        basic_stat=results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_annotated_tmp.bed",
        vcf_cols="header_annotated_vcf.txt",
        aa_af=expand(data_path + "AA_chr{chr}.afreq", chr=chromosomes),
        afr_af=expand(data_path + "AFR_chr{chr}.afreq", chr=chromosomes),
        eur_af=expand(data_path + "EUR_chr{chr}.afreq", chr=chromosomes),
        eas_af=expand(data_path + "EAS_chr{chr}.afreq", chr=chromosomes)
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_annotated.bed"
    params:
        archaic_genomes = " ".join(archaic_genomes)
    resources:
        mem_mb=64 * 1000,
        load=10

use rule aggregate_modern_variant_annotations_per_selected_segment as aggregate_modern_variant_annotations_per_known_desert with:
    input:
        modern_var=results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_annotated_tmp2.bed",
        basic_stat=results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_annotated_tmp.bed",
        vcf_cols="header_annotated_vcf.txt",
        aa_af=expand(data_path + "AA_chr{chr}.afreq", chr=chromosomes),
        afr_af=expand(data_path + "AFR_chr{chr}.afreq", chr=chromosomes),
        eur_af=expand(data_path + "EUR_chr{chr}.afreq", chr=chromosomes),
        eas_af=expand(data_path + "EAS_chr{chr}.afreq", chr=chromosomes)
    output:
        results_path + "ibdmix_{archaic_genome}/known_introgression_deserts_hg38_annotated.bed"
    params:
        archaic_genomes = " ".join(archaic_genomes)
    resources:
        mem_mb=64 * 1000,
        load=10

rule get_list_of_genes_overlapping_introgressed_segments:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/ibdmix_results_masked_denisovan_combined_" +
               str(int(minimum_length / 1000)) + "kb_" + str(lod_threshold) + "LOD.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz', '.bed.gz'),
        genes_to_ensembl=data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz', ''),
        genome_file=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        temp(results_path + "ibdmix_{archaic_genome}/background_list_of_genes_with_introgression.txt")
    params:
        bedtools=bedtools_path,
        aa_only = '-a',
        header='',
        slop=0
    shell:
        "scripts/get_gene_list_overlapping_introgressed_segments.sh -b {params.bedtools} -i {input.ibdmix} "
        "-g {input.genes} -e {input.genes_to_ensembl} -c {input.genome_file} -s {params.slop} -o {output} "
        "{params.aa_only} {params.header}"

rule get_list_of_genes_overlapping_pos_selected_introgressed_segments:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        genes_to_ensembl=data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz',''),
        genome_file=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        temp(results_path + "ibdmix_{archaic_genome}/foreground_list_of_genes_with_pos_selected_introgression.txt")
    params:
        bedtools=bedtools_path,
        aa_only= '',
        header='-h',
        slop=250000 #add 250kb to each site of an interval to get more genes
    shell:
        "scripts/get_gene_list_overlapping_introgressed_segments.sh -b {params.bedtools} "
        "-i <(awk -F '\\t' '{{if ($4 > $11) print $0}}' {input.ibdmix}) "
        "-g {input.genes} -e {input.genes_to_ensembl} -c {input.genome_file} -s {params.slop} -o {output} "
        "{params.aa_only} {params.header}"

rule get_list_of_genes_overlapping_neg_selected_introgressed_segments:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        genes_to_ensembl=data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz',''),
        genome_file=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        temp(results_path + "ibdmix_{archaic_genome}/foreground_list_of_genes_with_neg_selected_introgression.txt")
    params:
        bedtools=bedtools_path,
        aa_only= '',
        header='-h',
        slop=250000 #add 250kb to each site of an interval to get more genes
    shell:
        "scripts/get_gene_list_overlapping_introgressed_segments.sh -b {params.bedtools} "
        "-i <(awk -F '\\t' '{{if ($4 < $11) print $0}}' {input.ibdmix}) "
        "-g {input.genes} -e {input.genes_to_ensembl} -c {input.genome_file} -s {params.slop} -o {output} "
        "{params.aa_only} {params.header}"

use rule get_list_of_genes_overlapping_introgressed_segments as get_list_of_genes_overlapping_introgression_deserts with:
    input:
        ibdmix=results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        genes_to_ensembl=data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz',''),
        genome_file=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        temp(results_path + "ibdmix_{archaic_genome}/{superpopulation}_foreground_list_of_genes_in_deserts.txt")
    params:
        bedtools=bedtools_path,
        aa_only = '',
        header='',
        slop=0

use rule get_list_of_genes_overlapping_introgressed_segments as get_list_of_genes_overlapping_known_introgression_deserts with:
    input:
        ibdmix=results_path + "known_introgression_deserts_hg38_merged.bed",
        genes=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        genes_to_ensembl=data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz',''),
        genome_file=data_path + reference_path + 'genomefile_hg38.bed'
    output:
        temp(results_path + "foreground_list_of_genes_in_known_deserts.txt")
    params:
        bedtools=bedtools_path,
        aa_only = '',
        header='',
        slop=0

rule convert_gene_lists_bg:
    input:
        results_path + "ibdmix_{archaic_genome}/background_list_of_genes_with_introgression.txt"
    output:
        results_path + "ibdmix_{archaic_genome}/background_list_of_genes_with_introgression_converted.txt"
    conda:
        "../envs/biomart.yaml"
    shell:
        "scripts/look_up_gene_ids.R {input} {output}"

use rule convert_gene_lists_bg as convert_gene_lists_fg with:
    input:
        results_path+ "ibdmix_{archaic_genome}/foreground_list_of_genes_with_{segment_type}_selected_introgression.txt"
    output:
        results_path + "ibdmix_{archaic_genome}/foreground_list_of_genes_with_{segment_type}_selected_introgression_converted.txt"

use rule convert_gene_lists_bg as convert_gene_lists_deserts with:
    input:
        results_path + "ibdmix_{archaic_genome}/{superpopulation}_foreground_list_of_genes_in_deserts.txt"
    output:
        results_path + "ibdmix_{archaic_genome}/{superpopulation}_foreground_list_of_genes_in_deserts_converted.txt"

use rule convert_gene_lists_bg as convert_gene_lists_known_deserts with:
    input:
        results_path + "foreground_list_of_genes_in_known_deserts.txt"
    output:
        results_path + "foreground_list_of_genes_in_known_deserts_converted.txt"

rule get_gwas_hits_intersecting_putatively_selected_introgressed_regions:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_gwas_hits.bed"
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} intersect -a <(tail -n+2 {input.gwas}) -b {input.ibdmix} -wo > {output}"

use rule get_gwas_hits_intersecting_putatively_selected_introgressed_regions as get_gwas_hits_intersecting_introgression_deserts with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_gwas_hits.bed"
    params:
        bedtools = bedtools_path

use rule get_gwas_hits_intersecting_putatively_selected_introgressed_regions as get_gwas_hits_intersecting_introgression_desert_control_segments with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_gwas_hits.bed"
    params:
        bedtools = bedtools_path

use rule get_gwas_hits_intersecting_putatively_selected_introgressed_regions as get_gwas_hits_intersecting_control_regions with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_gwas_hits.bed"
    params:
        bedtools = bedtools_path
    wildcard_constraints:
        segment_type = "not_selected_control|non_introgressed_control"

rule get_independent_gwas_hits_per_putatively_selected_region:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_independent_gwas_hits.bed"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"
    shell:
        "scripts/clump_gwas_hits.R -c {input.gwas} -s {input.ibdmix} -a DISEASE.TRAIT -t {params.token} "
        "-o {output} -p {input.overlap}"

rule get_independent_gwas_hits_per_putatively_selected_region_parent_terms:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_independent_gwas_hits_parent_terms.bed"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"
    shell:
        "scripts/clump_gwas_hits.R -c {input.gwas} -s {input.ibdmix} -a Parent.term -t {params.token} "
        "-o {output} -p {input.overlap}"

use rule get_independent_gwas_hits_per_putatively_selected_region as get_independent_gwas_hits_per_desert with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_independent_gwas_hits.bed"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"

use rule get_independent_gwas_hits_per_putatively_selected_region_parent_terms as get_independent_gwas_hits_per_desert_parent_terms with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_independent_gwas_hits_parent_terms.bed"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"

use rule get_independent_gwas_hits_per_putatively_selected_region as get_independent_gwas_hits_per_desert_control_segments with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_independent_gwas_hits.bed"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"

use rule get_independent_gwas_hits_per_putatively_selected_region_parent_terms as get_independent_gwas_hits_per_desert_control_segments_parent_terms with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_independent_gwas_hits_parent_terms.bed"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"

use rule get_independent_gwas_hits_per_putatively_selected_region as get_independent_gwas_hits_per_control_regions with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_independent_gwas_hits.bed"
    wildcard_constraints:
        segment_type = "not_selected_control|non_introgressed_control"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"

use rule get_independent_gwas_hits_per_putatively_selected_region_parent_terms as get_independent_gwas_hits_per_control_regions_parent_terms with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
        gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.bed",
        overlap = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_gwas_hits.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_independent_gwas_hits_parent_terms.bed"
    wildcard_constraints:
        segment_type = "not_selected_control|non_introgressed_control"
    params:
        token=ldlink_token
    conda:
        "../envs/snpclip.yaml"

rule get_eqtls_overlapping_putatively_selected_regions:
    input:
        eqtls=data_path + reference_path + "GTEx_significant_independent_eqtls_all_tissues.bed",
        ibdmix= results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_iDAT_annotated.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_selected_neanderthal_segments_overlap_eQTLs.bed"
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} intersect -a {input.eqtls} -b {input.ibdmix} -wo > {output}"

use rule get_eqtls_overlapping_putatively_selected_regions as get_eqtls_overlapping_introgression_deserts with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_iDAT_annotated.bed",
        eqtls=data_path + reference_path + "GTEx_significant_independent_eqtls_all_tissues.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_novel_introgression_deserts_overlap_eQTLs.bed"
    params:
        bedtools = bedtools_path

use rule get_eqtls_overlapping_putatively_selected_regions as get_eqtls_overlapping_introgression_desert_control_segments with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_iDAT_annotated.bed",
        eqtls=data_path + reference_path + "GTEx_significant_independent_eqtls_all_tissues.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_introgression_deserts_new_control_segments_{n}_overlap_eQTLs.bed"
    params:
        bedtools = bedtools_path

use rule get_eqtls_overlapping_putatively_selected_regions as get_eqtls_overlapping_intersecting_control_regions with:
    input:
        ibdmix = results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_iDAT_{n}_annotated.bed",
        eqtls=data_path + reference_path + "GTEx_significant_independent_eqtls_all_tissues.bed"
    output:
        results_path + "ibdmix_{archaic_genome}/AMR_putatively_{segment_type}_neanderthal_segments_{n}_overlap_eQTLs.bed"
    params:
        bedtools = bedtools_path
    wildcard_constraints:
        segment_type = "not_selected_control|non_introgressed_control"
