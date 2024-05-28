"""
Extract ACB and ASW from 1KGP
"""
rule get_acb_asw_1kg:
    input:
        vcf=data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        tbi=data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
        panel=data_path + "1000G_phase3/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        vcf=data_path + 'target_individuals_chr{chr}.vcf.gz'
    threads: 8
    params:
        bcftools = bcftools_path
    shell:
        "{params.bcftools} view -S <(grep -E \"ACB|ASW\" {input.panel} | cut -f1) -Oz -o {output} --threads {threads} {input.vcf}"

rule index_acb_asw:
    input:
        data_path + 'target_individuals_chr{chr}.vcf.gz'
    output:
        data_path + "target_individuals_chr{chr}.vcf.gz.tbi"
    conda:
        "../envs/tabix.yaml"
    shell:
        "tabix {input}"

rule acb_asw_vcf_to_plink:
    input:
        vcf=data_path + 'target_individuals_chr{chr}.vcf.gz',
        fasta=data_path + reference_path + "hg38.fa"
    output:
        multiext(data_path + 'target_individuals_chr{chr}', ".bed", ".bim", ".fam")
    params:
        base = data_path + 'target_individuals_chr{chr}',
        maf = maf,
        geno = geno
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    threads: 8
    shell:
        "plink2 --vcf {input.vcf} --snps-only --chr 1-22 --allow-extra-chr --max-alleles 2 --memory {resources.mem_mb} "
        "--set-all-var-ids @:#\$1:\$2 --make-bed --out {params.base} --threads {threads} --ref-from-fa {input.fasta} "
        "--maf {params.maf} --geno {params.geno}"



