import gzip
import pandas as pd

rule create_rye_reference:
    output:
        data_path + merged_datasets_path + "pop2group.txt"
    shell:
        "echo -e 'Pop\tSubgroup\tGroup' > {output}; "
        "echo -e 'YRI\tWestAfrican\tWestAfrican' >> {output}; "
        "echo -e 'ESN\tWestAfrican\tWestAfrican' >> {output}; "
        "echo -e 'MSL\tWestAfrican\tWestAfrican' >> {output}; "
        "echo -e 'GWD\tWestAfrican\tWestAfrican' >> {output}; "
        "echo -e 'LWK\tEASTAfrican\tEastAfrican' >> {output}; "
        "echo -e 'CDX\tEastAsian\tEastAsian' >> {output}; "
        "echo -e 'CHB\tEastAsian\tEastAsian' >> {output}; "
        "echo -e 'CHS\tEastAsian\tEastAsian' >> {output}; "
        "echo -e 'KHV\tEastAsian\tEastAsian' >> {output}; "
        "echo -e 'JPT\tEastAsian\tEastAsian' >> {output}; "
        "echo -e 'Mayan\tMayan\tMayan' >> {output}; "
        "echo -e 'Pima\tPima\tPima' >> {output}; "
        "echo -e 'CEU\tWesternEuropean\tWesternEuropean' >> {output}; "
        "echo -e 'GBR\tWesternEuropean\tWesternEuropean' >> {output}; "
        "echo -e 'IBS\tIberian\tIberian' >> {output}; "
        "echo -e 'TSI\tSouthernEuropean\tSouthernEuropean' >> {output}; "
        "echo -e 'FIN\tNorthernEuropean\tNorthernEuropean' >> {output}; "

#
rule compute_pca_reference:
    input:
        multiext(data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned',
                 ".bed", ".fam", ".bim"),
    output:
        multiext(data_path + merged_datasets_path + 'reference_target_individuals_ALL_ld_pruned_ref',
                 ".eigenvec", ".eigenval", '.acount', ".eigenvec.allele")
    params:
        input_base = data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned',
        output_base = data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned_ref',
        pcs = M,
        maf = maf
    conda:
        "../envs/plink.yaml"
    threads: 16
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --bfile {params.input_base} --freq counts --pca {params.pcs} allele-wts --out {params.output_base} "
        "--threads {threads} --memory {resources.mem_mb} --maf 0.01 "
        "--keep <(cat {params.input_base}.fam | cut -f1,2 | grep -v -E '^0')"

# rule project_target_samples:
#     input:
#         target=multiext(data_path + 'target_individuals_chr{chr}', ".bim", ".bed", ".fam"),
#         acount=data_path + merged_datasets_path + "1KG_HGDP_ALL_shared_variants_w_target_ld_pruned.acount",
#         wts=data_path + merged_datasets_path + "1KG_HGDP_ALL_shared_variants_w_target_ld_pruned.eigenvec.allele"
#     output:
#         data_path + 'target_individuals_chr{chr}_projected.sscore'
#     conda:
#         "../envs/plink.yaml"
#     params:
#         target_prefix=data_path + 'target_individuals_chr{chr}',
#         out_prefix= data_path + 'target_individuals_chr{chr}' + "_projected",
#         cols=6 + M - 1
#     shell:
#         "plink2 --bfile {params.target_prefix} --read-freq {input.acount} --score {input.wts} 2 5 "
#         "header-read no-mean-imputation variance-standardize --score-col-nums 6-{params.cols} --out {params.out_prefix}"
# #
# rule aggregate_projections:
#     input:
#         expand(data_path + 'target_individuals_chr{chr}_projected.sscore', chr=chromosomes)
#     output:
#         data_path + merged_datasets_path + "target_samples_projected.sscore"
#     run:
#         scores = [pd.read_csv(score,sep='\t',header=0) for score in input]
#         allele_counts = scores[0].ALLELE_CT.values.copy()
#         for score in scores[1:]:
#             allele_counts += score.ALLELE_CT.values
#         projection = scores[0].copy()
#         projection.iloc[:, 4:] = projection.iloc[:, 4:].values * (projection.ALLELE_CT.values / allele_counts)[:, np.newaxis]
#         projection.sort_values('IID',inplace=True)
#         for score in scores[1:]:
#             score.sort_values('IID',inplace=True)
#             projection.iloc[:, 4:] += score.iloc[:, 4:].values * (score.ALLELE_CT.values / allele_counts)[:,
#                                                                       np.newaxis]
#             projection.ALLELE_CT += score.ALLELE_CT.values
#             projection.NAMED_ALLELE_DOSAGE_SUM += score.NAMED_ALLELE_DOSAGE_SUM.values
#         projection.to_csv(output[0], header=True, sep='\t', index=False)
# #
# #
rule project_samples:
    input:
        samples=multiext(data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned", ".bim", ".bed", ".fam"),
        acount=data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned_ref.acount",
        wts=data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned_ref.eigenvec.allele"
    output:
        data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned.sscore"
    conda:
        "../envs/plink.yaml"
    params:
        prefix=data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned",
        cols=6 + M - 1
    shell:
        "plink2 --bfile {params.prefix} --read-freq {input.acount} --score {input.wts} 2 5 "
        "header-read no-mean-imputation variance-standardize --score-col-nums 6-{params.cols} --out {params.prefix}"

rule format_eigenvector_file_for_rye:
    input:
        data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned.sscore"
    output:
        data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned.eigenvec"
    shell:
        "cat {input} | sed 's/_AVG//g' | cut -f1,2,5- > {output}"
#
rule copy_eigenval:
    input:
        eigenval = data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned_ref.eigenval"
    output:
        data_path + merged_datasets_path  + "reference_target_individuals_ALL_ld_pruned.eigenval"
    shell:
        "cp {input} {output}"

# global ancestry proportions using rye
rule determine_ancestry_proportions:
    input:
        eigenvec=data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned.eigenvec',
        eigenval=data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned.eigenval',
        pop2group=data_path + merged_datasets_path + "pop2group.txt"
    output:
        f"{data_path}{merged_datasets_path}ancestry_proportions-{M}.{K}.Q",
        f"{data_path}{merged_datasets_path}ancestry_proportions-{M}.fam"
    params:
        rye = rye_path,
        out = data_path + merged_datasets_path + "ancestry_proportions",
        pcs = M
    threads: 16
    conda:
        "../envs/rye.yaml"
    shell:
        "{params.rye} --eigenval={input.eigenval} --eigenvec={input.eigenvec} --pop2group={input.pop2group} "
        "--output={params.out} --pcs={params.pcs} --threads={threads}"
#
rule select_african_american_individuals:
    input:
        qfile=f"{data_path}{merged_datasets_path}ancestry_proportions-{M}.{K}.Q",
        fam=f"{data_path}{merged_datasets_path}ancestry_proportions-{M}.fam",
        ref=data_path + "1000G_phase3/reference_panel_filtered.txt"
    output:
        data_path + "AA_sample_ids.txt"
    params:
        min_afr = min_afr_component,
        max_afr = max_afr_component,
        min_eur = min_eur_component,
        max_eur = max_eur_component,
        min_sum_afr_eur = min_afr_eur_component_combined
    shell:
        "awk -F '\\t' '{{if ($2 + $3 >= {params.min_afr} && $7 + $8 + $9 + $10 >= {params.min_eur} && "
        "$2 + $3 + $7 + $8 + $9 + $10 >= {params.min_sum_afr_eur}) print $1}}' {input.qfile} |"
        " grep -v -f <(cut -f1 {input.ref}) > {output}"

# rule select_individuals_from_1kgp_and_merge_with_aa_individuals:
#     input:
#         ref_individuals=data_path + "1000G_phase3/reference_panel.txt",
#         aa_individuals=data_path + "AA_sample_ids.txt",
#         fam= data_path + merged_datasets_path + "1KG_HGDP_target_chr{chr}.fam"
#     output:
#         temp("african_american_and_ref_individuals_chr{chr}.txt")
#     shell:
#         "grep -f <(cut -f1 {input.ref_individuals}) {input.fam} | cut -d ' ' -f1,2 > {output}; "
#         "grep -f {input.aa_individuals} {input.fam} | cut -d ' ' -f1,2 >> {output}"
#
# rule extract_african_american_individuals_and_1kgp:
#     input:
#         aa_1kgp_indv = "african_american_and_ref_individuals_chr{chr}.txt",
#         plink=multiext(data_path + merged_datasets_path + "1KG_HGDP_target_chr{chr}", ".bed", ".bim", ".fam")
#     output:
#         temp(multiext(data_path + "african_american_and_ref_individuals_chr{chr}", ".bed", ".bim", ".fam"))
#     conda:
#         "../envs/plink.yaml"
#     threads: 16
#     resources:
#         mem_mb= lambda wildcards,threads: threads * mem_gb_cpu * 1000
#     params:
#         input_base = data_path + merged_datasets_path + "1KG_HGDP_target_chr{chr}",
#         output_base = data_path + "african_american_and_ref_individuals_chr{chr}"
#     shell:
#         "plink2 --bfile {params.input_base} --keep {input.aa_1kgp_indv} --make-bed --out {params.output_base} "
#         "--threads {threads} --memory {resources.mem_mb}"
#
# # convert plink files to vcf
# rule plink_to_vcf_target:
#     input:
#         fasta=data_path + reference_path + "hg38.fa",
#         plink=multiext(data_path + "african_american_and_ref_individuals_chr{chr}", ".bed", ".bim", ".fam")
#     output:
#         temp(data_path + "african_american_and_ref_individuals_chr{chr}.vcf.gz")
#     params:
#         base = data_path + "african_american_and_ref_individuals_chr{chr}"
#     conda:
#         "../envs/plink.yaml"
#     resources:
#         mem_mb= lambda wildcards,threads: threads * mem_gb_cpu * 1000
#     threads: 8
#     shell:
#         "plink2 --bfile {params.base} --export vcf bgz id-paste=iid --out {params.base} --threads {threads} "
#         "--memory {resources.mem_mb} --ref-from-fa {input.fasta}"

# train flare model for local ancestry inference on chromsome 1
rule train_flare_model:
    input:
        ref_vcf = data_path + "1000G_phase3/REF.chr1.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        ref_panel = data_path + "1000G_phase3/reference_panel_filtered.txt",
        map = data_path + genetic_map_path + "plink.formatted.chr1.GRCh38.map",
        target_vcf=data_path + 'target_individuals_chr1.vcf.gz',
        aa_individuals = data_path + "AA_sample_ids.txt",
    output:
        multiext(results_path + flare_output + "african_american_and_ref_individuals_chr1",
            ".model", ".anc.vcf.gz")
    params:
        base = results_path + flare_output + "african_american_and_ref_individuals_chr1",
        flare = flare_path,
        seed = 42,
        mem_gb= lambda wildcards,threads: threads * mem_gb_cpu
    threads: 16
    resources:
        mem_mb= lambda wildcards,threads: threads * mem_gb_cpu * 1000
    shell:
        "java -Xmx{params.mem_gb}g -jar {params.flare} ref={input.ref_vcf} "
        "ref-panel={input.ref_panel} gt={input.target_vcf} gt-samples={input.aa_individuals} map={input.map} "
        "out={params.base} nthreads={threads} seed={params.seed}"

# paint chromosomes using previously trained flare model
rule run_flare:
    input:
        ref_vcf = data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        ref_panel = data_path + "1000G_phase3/reference_panel_filtered.txt",
        map = data_path + genetic_map_path + "plink.formatted.chr{chr}.GRCh38.map",
        model = results_path + flare_output + "african_american_and_ref_individuals_chr1.model",
        target_vcf=data_path + 'target_individuals_chr{chr}.vcf.gz',
        aa_individuals = data_path + "AA_sample_ids.txt",
    output:
        results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc.vcf.gz"
    wildcard_constraints:
        chr="[1-9][0-9]?|1[0-9]"
    params:
        base = results_path + flare_output +"african_american_and_ref_individuals_chr{chr}",
        flare = flare_path,
        seed = 42,
        mem_gb = lambda wildcards,threads: threads * mem_gb_cpu
    threads: 16
    resources:
        mem_mb= lambda wildcards,threads: threads * mem_gb_cpu * 1000,
        load=10
    shell:
        "java -Xmx{params.mem_gb}g -jar {params.flare} ref={input.ref_vcf} "
        "ref-panel={input.ref_panel} gt={input.target_vcf} gt-samples={input.aa_individuals} map={input.map} "
        "out={params.base} nthreads={threads} seed={params.seed} em=false model={input.model}"


def get_params_extract_local_ancestry(wildcards):
    if wildcards.phase == '0':
        params = {"include0": "'GT=\\\"0|1\\\" | GT=\\\"0|0\\\"'",
                  "include1": "'GT=\\\"1|0\\\" | GT=\\\"1|1\\\"'",
                  "fields0": "'[%CHROM\\t%POS\\t%REF\\t%AN1\\t%SAMPLE\\n]'",
                  "fields1": "'[%CHROM\\t%POS\\t%ALT\\t%AN1\\t%SAMPLE\\n]'"}
    elif wildcards.phase == '1':
        params = {"include0": "'GT=\\\"1|0\\\" | GT=\\\"0|0\\\"'",
                  "include1": "'GT=\\\"0|1\\\" | GT=\\\"1|1\\\"'",
                  "fields0": "'[%CHROM\\t%POS\\t%REF\\t%AN1\\t%SAMPLE\\n]'",
                  "fields1": "'[%CHROM\\t%POS\\t%ALT\\t%AN1\\t%SAMPLE\\n]'"}
    return params

rule extract_local_ancestry:
    input:
        lai=results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc.vcf.gz",
    output:
        phase0=results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase0.{archaic_genome}.bed",
        phase1=results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase1.{archaic_genome}.bed"
    shell:
        "scripts/extract_local_ancestry_information.py --lai {input.lai} --output_phase0 {output.phase0} "
        "--output_phase1 {output.phase1}"


rule compute_local_ancestry_coverage_per_window_and_phase:
    input:
        genomefile=data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed",
        lai=results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_pos.phase{phase}.{archaic_genome}.bed",
    output:
        temp(results_path + flare_output + \
        "african_american_and_ref_individuals_chr{chr}.anc_per_window{windowsize}_s_{stepsize}.phase{phase}.{archaic_genome}.bed")
    params:
        bedtools = bedtools_path,
    wildcard_constraints:
        phase='0|1',
        archaic_genome = '|'.join(neanderthal_genomes)
    threads: 16
    shell:
        "{params.bedtools} intersect -a {input.genomefile} -b {input.lai} -wo | sort -k1,1 -k2,2n -k3,3n -k7 "
        "--parallel={threads} | {params.bedtools} groupby -g 1,2,3,7 -c 8 -o count > {output}"

rule compute_local_ancestry_frequencies_per_window:
    input:
        lai_0 = results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_window{windowsize}_s_{stepsize}.phase0.{archaic_genome}.bed",
        lai_1= results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_window{windowsize}_s_{stepsize}.phase1.{archaic_genome}.bed"
    output:
        results_path + flare_output + "african_american_and_ref_individuals_chr{chr}.anc_per_window{windowsize}_s_{stepsize}.{archaic_genome}.bed"
    wildcard_constraints:
        windowsize="[0-9]+",
        stepsize="[0-9]+",
        archaic_genome='|'.join(neanderthal_genomes)
    params:
        bedtools = bedtools_path
    shell:
        "cat {input.lai_0} {input.lai_1} | sort -k1,1 -k2,2n -k3,3n -k4 | "
        "{params.bedtools} groupby -g 1-4 -c 5 -o sum  > {output}"
