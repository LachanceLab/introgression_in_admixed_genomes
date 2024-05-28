"""
Phase All of Us genotype data
"""

rule download_beagle:
    output:
        beagle_url.split('/')[-1]
    params:
        url = beagle_url
    shell:
        "wget -q {params.url}"

rule phase_all_of_us:
    input:
        vcf=GS.remote(GS_PREFIX + '/' + all_of_us_vcf),
        tbi=GS.remote(GS_PREFIX + '/' + all_of_us_vcf + ".tbi"),
        map=data_path + genetic_map_path + "plink.formatted.chr{chr}.GRCh38.map",
        beagle=beagle_url.split('/')[-1],
        samples_to_exclude='non_target_samples_all_of_us.txt'
    output:
        phased=data_path + 'target_individuals_chr{chr}.vcf.gz',
        gz=temp(GS_PREFIX + '/' + all_of_us_vcf.replace('.bgz', '.gz')),
        tbi=temp(GS_PREFIX + '/' + all_of_us_vcf.replace('.bgz', '.gz') + ".tbi")
    conda:
        "../envs/beagle.yaml"
    params:
        jar = beagle_url.split('/')[-1],
        output_prefix=data_path + 'target_individuals_chr{chr}',
        vcfgz=GS_PREFIX + '/' + all_of_us_vcf.replace('.bgz', '.gz'),
        tbigz=GS_PREFIX + '/' + all_of_us_vcf.replace('.bgz', '.gz') + ".tbi"
    threads: 96
    resources:
        mem_gb = lambda wildcards, threads: threads * mem_gb_cpu
    benchmark:
        "benchmark_beagle_all_of_us_chr{chr}.txt"
    shell:
        "mv {input.vcf} {params.vcfgz}; mv {input.tbi} {params.tbigz}; "
        "java -Xmx{resources.mem_gb}g -jar {params.jar} nthreads={threads} gt={params.vcfgz} "
        "map={input.map} out={params.output_prefix} excludesamples={input.samples_to_exclude}"

rule index_phased_vcf_all_of_us:
    input:
        data_path + 'target_individuals_chr{chr}.vcf.gz'
    output:
        data_path + 'target_individuals_chr{chr}.vcf.gz.tbi'
    conda:
        "../envs/tabix.yaml"
    shell:
        "tabix {input}"

rule target_vcf_to_plink:
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
