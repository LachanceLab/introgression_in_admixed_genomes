import gzip
import pandas as pd
import numpy as np

rule download_hg19_to_hg38_chain:
    output:
        data_path + reference_path + "hg19ToHg38.over.chain.gz"
    params:
        data_path = data_path + reference_path,
        url = hg19_to_hg38_chain_url
    shell:
        "wget -q -P {params.data_path} {params.url}"

# Download 1KGP Phase3 VCF files and sample-population mapping file. Reference genome is hg19
rule download_1kgp_phase3:
    output:
        vcf = temp(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_hg19.vcf.gz"),
        tbi = temp(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_hg19.vcf.gz.tbi")
    params:
        vcf = phase3_1KG_base_url + "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        tbi = phase3_1KG_base_url + "ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
    shell:
        "wget -q -O {output.vcf} {params.vcf}; " 
        "wget -q -O {output.tbi} {params.tbi};"


rule lift_1kgp_from_hg19_to_hg38:
    input:
        vcf = data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_hg19.vcf.gz",
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz",
        reference=data_path + reference_path + hg38_fasta_url.split('/')[-1].replace('.gz', '')
    output:
        vcf = temp(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"),
        unmapped = temp(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.unmap")
    params:
        bcftools = bcftools_path
    conda:
        '../envs/crossmap.yaml'
    shell:
        "CrossMap.py vcf --chromid l {input.chain} <({params.bcftools} view -v snps {input.vcf}) {input.reference} {output.vcf}"

rule compress_and_index_1kgp:
    input:
        data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf"
    output:
        vcf=temp(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_tmp.vcf.gz"),
        tbi=temp(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_tmp.vcf.gz.tbi")
    conda:
        "../envs/tabix.yaml"
    params:
        tmp_dir = data_path + "1000G_phase3/",
        bcftools = bcftools_path
    shell:
        "{params.bcftools} sort -T {params.tmp_dir} -Oz -o {output.vcf} {input}; tabix {output.vcf}"

rule extract_chromosome_vcf_1kgp:
    input:
        vcf=expand(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_tmp.vcf.gz",
                   chr=chromosomes),
        tbi=expand(data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes_tmp.vcf.gz.tbi",
                   chr=chromosomes)
    output:
        vcf=data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        tbi=data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
    params:
        chrom = '{chr}',
        tmp_dir= data_path + "1000G_phase3/",
        bcftools = bcftools_path
    conda:
        "../envs/tabix.yaml"
    shell:
        "{params.bcftools} concat -a -r 'chr{params.chrom}' {input.vcf} | {params.bcftools} sort -T {params.tmp_dir} -Oz -o {output.vcf}; "
        "tabix {output.vcf}"

# get file that maps sample IDs to populations
rule download_1kgp_panel_file:
    output:
        samples = data_path + "1000G_phase3/integrated_call_samples_v3.20130502.ALL.panel"
    params:
        samples = phase3_1KG_base_url + "integrated_call_samples_v3.20130502.ALL.panel",
        data_path= data_path
    shell:
        "mkdir -p {params.data_path}1000G_phase3; "
        "wget -q -O {output.samples} {params.samples}"

# download accessible genome mask 1000G --> positions should be included --> hg19 --> liftover to hg38
rule download_accessible_genome_mask:
    output:
        temp(data_path + mask_path + 'tmp_' + accessible_genome_mask_url.split('/')[-1])
    params:
        url = accessible_genome_mask_url
    shell:
        "wget -qO- {params.url} | cut -f1-3 | sed 's/^chr//' > {output}"

rule lift_1kgp_mask_from_hg19_to_hg38:
    input:
        bed=data_path + mask_path + 'tmp_' + accessible_genome_mask_url.split('/')[-1],
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz"
    output:
        data_path + mask_path + accessible_genome_mask_url.split('/')[-1]
    conda:
        "../envs/crossmap.yaml"
    shell:
        'CrossMap.py bed --chromid l {input.chain} {input.bed} {output}'

rule download_hgdp_wgs:
    output:
        vcf=temp(data_path + "hgdp/hgdp_wgs.20190516.full.chr{chr}.vcf.gz"),
        tbi=temp(data_path + "hgdp/hgdp_wgs.20190516.full.chr{chr}.vcf.gz.tbi")
    params:
        url_vcf=hgdp_wgs_url,
        url_tbi=hgdp_wgs_url + '.tbi'
    retries: 2
    shell:
        "wget -q -O {output.vcf} {params.url_vcf}; "
        "wget -q -O {output.tbi} {params.url_tbi}"

def get_urls_download_archaic_genome(wildcards):
    urls = dict()
    if wildcards.archaic_genome == 'Chagyrskaya':
        urls['vcf'] = archaic_genomes_mpg_base_url_chagyrskaya + f"VCF/chr{wildcards.chr}.noRB.vcf.gz"
        urls['tbi'] = archaic_genomes_mpg_base_url_chagyrskaya + f"VCF/chr{wildcards.chr}.noRB.vcf.gz.tbi"
        urls['filter_bed'] = archaic_genomes_mpg_base_url_chagyrskaya + f"FilterBed/chr{wildcards.chr}_mask.bed.gz"
    elif wildcards.archaic_genome == 'AltaiNeandertal':
        urls['vcf'] = archaic_genomes_mpg_base_url_altai_old + f"{wildcards.archaic_genome}/VCF/AltaiNea.hg19_1000g.{wildcards.chr}.mod.vcf.gz"
        urls['filter_bed'] = "https://bioinf.eva.mpg.de/altai_minimal_filters/AltaiNea.map35_50.MQ30.Cov.indels.TRF.bed.bgz"
    else:
        urls['vcf'] = archaic_genomes_mpg_base_url + f"VCF/{wildcards.archaic_genome}/chr{wildcards.chr}_mq25_mapab100.vcf.gz"
        urls['tbi'] = archaic_genomes_mpg_base_url + f"VCF/{wildcards.archaic_genome}/chr{wildcards.chr}_mq25_mapab100.vcf.gz.tbi"
        urls['filter_bed'] = archaic_genomes_mpg_base_url + f"FilterBed/{wildcards.archaic_genome}/chr{wildcards.chr}_mask.bed.gz"
    return urls

# Download archaic genomes (VCFs) and bed files with regions which to include. Reference genome is hg19
rule download_archaic_genome:
    output:
        vcf=temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_hg19.vcf.gz"),
        tbi=temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_hg19.vcf.gz.tbi"),
        filter_bed=temp(data_path + "archaic_genomes/{archaic_genome}/filter_bed/chr{chr}_mask_hg19.bed.gz")
    params:
        urls = get_urls_download_archaic_genome,
        data_path= data_path
    conda:
        '../envs/tabix.yaml'
    wildcard_constraints:
        archaic_genome="Altai|Vindija33.19|Chagyrskaya|Denisova"
    shell:
        "mkdir -p {params.data_path}archaic_genomes/{wildcards.archaic_genome}/vcf/; "
        "mkdir -p {params.data_path}archaic_genomes/{wildcards.archaic_genome}/filter_bed/; "
        "wget -q -O {output.vcf} {params[urls][vcf]}; "
        "wget -q -O {output.tbi} {params[urls][tbi]}; "
        "wget -q -O {output.filter_bed} {params[urls][filter_bed]}; "

rule download_archaic_genome_altai_old:
    output:
        vcf=temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_hg19.vcf.gz"),
        tbi=temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_hg19.vcf.gz.tbi"),
        filter_bed=temp(data_path + "archaic_genomes/{archaic_genome}/filter_bed/chr{chr}_mask_hg19.bed.gz")
    params:
        urls = get_urls_download_archaic_genome,
        data_path= data_path
    conda:
        '../envs/tabix.yaml'
    wildcard_constraints:
        archaic_genome='AltaiNeandertal'
    shell:
        "mkdir -p {params.data_path}archaic_genomes/{wildcards.archaic_genome}/vcf/; "
        "mkdir -p {params.data_path}archaic_genomes/{wildcards.archaic_genome}/filter_bed/; "
        "wget -q -O {output.vcf} {params[urls][vcf]}; "
        "tabix {output.vcf}; "
        "wget -q -O - {params[urls][filter_bed]} | zegrep -w \"^{wildcards.chr}\" | gzip > {output.filter_bed}; "

use rule lift_1kgp_from_hg19_to_hg38 as lift_archaic_genomes_from_hg19_to_hg38 with:
    input:
        vcf = data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_hg19.vcf.gz",
        tbi = data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_hg19.vcf.gz.tbi",
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz",
        reference=data_path + reference_path + hg38_fasta_url.split('/')[-1].replace('.gz', '')
    output:
        vcf = temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf"),
        unmapped = temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.unmap")

rule compress_and_index_archaic_genomes:
    input:
        data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf"
    output:
        vcf=temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_tmp.vcf.gz"),
        tbi=temp(data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100_tmp.vcf.gz.tbi")
    conda:
        "../envs/tabix.yaml"
    params:
        tmp_dir = data_path + "archaic_genomes/{archaic_genome}/",
        bcftools= bcftools_path
    shell:
        "{params.bcftools} sort -T {params.tmp_dir} -Oz -o {output.vcf} {input}; tabix {output.vcf}"

def get_lifted_over_archaic_chromosomes(wildcards):
    return [data_path + f"archaic_genomes/{wildcards.archaic_genome}/vcf/chr{chrom}_mq25_mapab100_tmp.vcf.gz" for chrom in chromosomes]

def get_lifted_over_archaic_chromosomes_indices(wildcards):
    return [data_path + f"archaic_genomes/{wildcards.archaic_genome}/vcf/chr{chrom}_mq25_mapab100_tmp.vcf.gz.tbi" for chrom in chromosomes]

rule extract_chromosome_vcf_archaic_genome:
    input:
        vcf=get_lifted_over_archaic_chromosomes,
        tbi=get_lifted_over_archaic_chromosomes_indices
    output:
        vcf=data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz",
        tbi=data_path + "archaic_genomes/{archaic_genome}/vcf/chr{chr}_mq25_mapab100.vcf.gz.tbi"
    params:
        chrom = '{chr}',
        tmp_dir= data_path + "archaic_genomes/{archaic_genome}/",
        bcftools= bcftools_path
    conda:
        "../envs/tabix.yaml"
    shell:
        "{params.bcftools} concat -a -r 'chr{params.chrom}' {input.vcf} | {params.bcftools} sort -T {params.tmp_dir} -Oz -o {output.vcf}; tabix {output.vcf} "

use rule lift_1kgp_mask_from_hg19_to_hg38 as lift_archaic_genomes_mask_from_hg19_to_hg38 with:
    input:
        bed=data_path + "archaic_genomes/{archaic_genome}/filter_bed/chr{chr}_mask_hg19.bed.gz",
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz"
    output:
        temp(data_path + "archaic_genomes/{archaic_genome}/filter_bed/chr{chr}_mask.bed")

def get_masks(wildcards):
    return [data_path + f"archaic_genomes/{wildcards.archaic_genome}/filter_bed/chr{chrom}_mask.bed" for chrom in chromosomes]

rule compress_archaic_mask:
    input:
        get_masks
    output:
        data_path + "archaic_genomes/{archaic_genome}/filter_bed/chr{chr}_mask.bed.gz"
    conda:
        "../envs/tabix.yaml"
    params:
        chrom = '{chr}'
    shell:
        "grep \"chr{params.chrom}\" {input} | cut -d ':' -f2 | sort -k1,1 -k2,2n  | bgzip > {output}"

rule download_ensemble_uniprot_table:
    output:
        temp(data_path + reference_path + ensembl_uniprot_url.split('/')[-1])
    params:
        url=ensembl_uniprot_url,
        output_dir=data_path + reference_path
    shell:
        "wget -q -P {params.output_dir} {params.url}"

rule decompress_ensemble_uniprot_table:
    input:
        data_path + reference_path + ensembl_uniprot_url.split('/')[-1]
    output:
        data_path + reference_path + ensembl_uniprot_url.split('/')[-1].replace('.gz', '')
    shell:
        "gunzip {input}"

rule download_string:
    output:
        temp(data_path + reference_path + string_url.split('/')[-1])
    params:
        url = string_url,
        output_dir = data_path + reference_path
    shell:
        "wget -q -P {params.output_dir} {params.url}"

rule decompress_string:
    input:
        data_path + reference_path + string_url.split('/')[-1]
    output:
        data_path + reference_path + string_url.split('/')[-1].replace('.gz', '')
    shell:
        "gunzip -c {input} | tail -n+2 | sed 's/9606.//g' > {output}"

rule download_gwas_efo_trait_mapping:
    output:
        data_path + reference_path + 'gwas_catalog_trait_mapping.tab'
    params:
        url = gwas_efo_trait_mapping_url
    shell:
        "wget -q -O {output} {params.url}"

rule download_gwas_catalog:
    output:
        temp(data_path + reference_path + "gwas_catalog_assocations_hg38.tab.tmp")
    params:
        url = gwas_catalog_url
    shell:
        # select only genome-wide significant hits
        "wget -q -O {output} {params.url}"

rule get_header_gwas_catalog:
    input:
        data_path + reference_path + "gwas_catalog_assocations_hg38.tab.tmp"
    output:
        temp(data_path + reference_path + "header.gwas")
    shell:
        "head -n1 {input} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $12, $13, $13\"_1\", $8, $14, $15, $16, $17, "
        "$18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $31, $32, $33, $34, $35, $36, $37, $7, $8, $38}}'> {output}"

rule filter_gwas_hits:
    input:
        header=data_path + reference_path + "header.gwas",
        gwas=data_path + reference_path + "gwas_catalog_assocations_hg38.tab.tmp"
    output:
        data_path + reference_path + "gwas_catalog_assocations_hg38.bed.tmp2"
    shell:
        "cat {input.header} > {output}; "
        "awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($28 < 5e-8 && $13 !~ /;/ && $12 !~ /x/) print \"chr\"$12, $13, $13 +1, $8, $14, $15, $16, $17, "
        "$18, $19, $20, $21, $22, $23, $24, $25, $26, $27, $28, $29, $31, $32, $33, $34, $35, $36, $37, $7, $8, $38}}' "
        "{input.gwas} | grep -E \"^chr[0-9]\" | sort -k1,1 -k2,2n >> {output}"

rule add_gwas_parent_terms:
    input:
        data_path + reference_path + "gwas_catalog_assocations_hg38.bed.tmp2",
        data_path+ reference_path + 'gwas_catalog_trait_mapping.tab'
    output:
        data_path + reference_path + "gwas_catalog_assocations_hg38.bed"
    run:
        gwas_catalog = pd.read_csv(input[0], header=0, sep='\t')
        gwas_trait_mapping = pd.read_csv(input[1],header=0,sep='\t')
        gwas_trait_mapping.loc[:, ['Disease trait', 'Parent term']].drop_duplicates(inplace=True)
        gwas_catalog_extended = gwas_catalog.set_index('DISEASE/TRAIT').join(
            gwas_trait_mapping.loc[:, ['Disease trait', 'Parent term']].set_index('Disease trait'))
        gwas_catalog_extended.reset_index(inplace=True)
        columns = gwas_catalog_extended.columns.values.tolist()
        columns_reordered = columns[1:]
        columns_reordered.append(columns[0])
        gwas_catalog_extended = gwas_catalog_extended.loc[:, columns_reordered]
        gwas_catalog_extended.drop_duplicates(["SNPS", "STUDY ACCESSION", "DISEASE/TRAIT", "Parent term"],
            inplace=True)
        gwas_catalog_extended.to_csv(output[0], header=True, index=False, sep='\t')

# rule get_dbSNP:
#     output:
#         temp(data_path + reference_path + dbSNP_url.split('/')[-1])
#     params:
#         url = dbSNP_url,
#         base_path = data_path + reference_path
#     shell:
#         "wget -q -P {params.base_path} {params.url}"
#
# rule convert_dbSNP_to_bed:
#     input:
#         data_path + reference_path + dbSNP_url.split('/')[-1]
#     output:
#         temp(data_path + reference_path + dbSNP_url.split('/')[-1].replace('.bb', '.bed'))
#     conda:
#         '../envs/bigbedtobed.yaml'
#     shell:
#         "bigBedToBed {input} {output}"


# rule update_coordinates_gwas_catalog_to_hg19:
#     input:
#         gwas = data_path + reference_path + "gwas_catalog_assocations_hg38.tab",
#         dbsnp_hg19 =  data_path + reference_path + dbSNP_url.split('/')[-1].replace('.bb', '.bed'),
#     output:
#         data_path + reference_path + "gwas_catalog_assocations_hg19.bed"
#     shell:
#         'scripts/update_coords_gwas_catalog_to_hg19.py --gwas {input.gwas} --dbsnp_hg19 {input.dbsnp_hg19} -o {output}'

# download map of recent segmental duplications in human genome
rule download_segmental_duplication_map:
    output:
        data_path + mask_path + segmental_duplications_url.split('/')[-1].split('.txt.gz')[0] + '.bed'
    params:
        download_filename = segmental_duplications_url.split('/')[-1],
        url = segmental_duplications_url,
        data_path = data_path,
        mask_path = mask_path,
        bedtools = bedtools_path
    shell:
        "wget -q -P {params.data_path}{params.mask_path} {params.url}; "
        "zcat {params.data_path}{params.mask_path}{params.download_filename} | cut -f2-4 | "
        "grep -P '^chr[0-9][0-9]?\t' | {params.bedtools} merge -i - | "
        "sort -k1,1 -k2,2n > {output}; rm {params.data_path}{params.mask_path}{params.download_filename}"

# download simple repeats identified by TRF
rule download_simple_repeat_map:
    output:
        data_path + mask_path + simple_repeat_map_url.split('/')[-1].split('.gz')[0]
    params:
        url = simple_repeat_map_url,
        data_path= data_path,
        mask_path = mask_path,
        bedtools = bedtools_path
    shell:
        "wget -q -P {params.data_path}{params.mask_path} {params.url}; "
        "zcat {output}.gz | cut -f2-4 | grep -P '^chr[0-9][0-9]?\t' | {params.bedtools} merge -i - | "
        "sort -k1,1 -k2,2n > {output}; rm {output}.gz"

rule download_gaps:
    output:
        data_path + mask_path + reference_gap_url.split('/')[-1].split('.gz')[0]
    params:
        url=reference_gap_url,
        data_path=data_path,
        mask_path=mask_path,
        bedtools= bedtools_path
    shell:
        "wget -q -O- {params.url} | gunzip | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $2, $3, $4}}' | "
        "grep -P -w 'chr[0-9][0-9]?' | sort -k1,1 -k2,2 > {output}"

# chrom txStart txEnd gene_id strand cdStart cdsEnd exonCount exonStarts exonEnds
rule download_annotated_genes:
    output:
        bed=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz', '.bed.gz'),
        tbi=data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz', '.bed.gz.tbi')
    params:
        url = knownGene_url,
        data_path = data_path,
        reference_path = reference_path,
        tmp_file= data_path + reference_path + knownGene_url.split('/')[-1]
    conda:
        "../envs/tabix.yaml"
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}; zcat {params.tmp_file} | "
        "awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $2, $4, $5, $1, $3, $6, $7, $8, $9, $10}}'| "
        "grep -E -w '^chr[0-9][0-9]?' | sort -k1,1 -k2,2n | bgzip > {output.bed}; "
        "tabix {output.bed}; rm {params.tmp_file}"

rule ucsc_gene_id_to_ensembl_mapping:
    output:
        data_path + reference_path + knownToEnsembl_url.split('/')[-1].replace('.gz', '')
    params:
        url=knownToEnsembl_url,
        output_dir=data_path + reference_path
    shell:
        "wget -q -P {params.output_dir} {params.url}; gunzip {output}.gz"

# chrom cdsstart cdsend cdsid
rule extract_CDS_from_knownGenes_table:
    input:
        data_path + reference_path + knownGene_url.split('/')[-1].replace('.txt.gz', '.bed.gz')
    output:
        data_path + reference_path + 'knownCDS.bed'
    run:
        cds = open(output[0], 'w')
        with gzip.open(input[0],'r') as kg:
            for line in kg:
                chrom, start, end, gene, strand, _, _, n_exons, exon_starts, exon_ends = line.strip().split(b'\t')
                i = 1
                for cds_start, cds_end in zip(exon_starts.split(b',')[:-1],exon_ends.split(b',')[:-1]):
                    cds.write(f'{str(chrom)[2:-1]}\t{str(cds_start)[2:-1]}\t{str(cds_end)[2:-1]}\t{str(gene)[2:-1]}.cds_' +
                              f'{i}\n')
                    i += 1
        kg.close()
        cds.close()

rule index_cds_bed:
    input:
        data_path + reference_path + 'knownCDS.bed'
    output:
        bed=data_path + reference_path + 'knownCDS.bed.gz',
        tbi=data_path + reference_path + 'knownCDS.bed.gz.tbi'
    conda:
        '../envs/tabix.yaml'
    shell:
        "sort -k1,1 -k2,2n {input} > {input}.tmp; mv {input}.tmp {input}; bgzip {input}; tabix {output.bed}"


# chrom start end id gerp
rule download_gerp_scores:
    output:
        temp(data_path + reference_path + gerp_url.split('/')[-1].replace('.bw', '.bed'))
    params:
        url = gerp_url,
        data_path = data_path,
        reference_path = reference_path,
        prefix = gerp_url.split('/')[-1].split(".bw")[0]
    conda:
        "../envs/bigwigtowig.yaml"
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}; "
        "bigWigToWig {params.data_path}{params.reference_path}{params.prefix}.bw "
        "{params.data_path}{params.reference_path}{params.prefix}.wig; "
        "wig2bed < {params.data_path}{params.reference_path}{params.prefix}.wig | "
        "grep -E -w '^chr[0-9][0-9]?' | sort -k1,1 -k2,2n > {output}; rm {params.data_path}{params.reference_path}{params.prefix}.bw; "
        "rm {params.data_path}{params.reference_path}{params.prefix}.wig"

rule compress_and_index_gerp_scores:
    input:
        data_path + reference_path + gerp_url.split('/')[-1].replace('.bw','.bed')
    output:
        bed=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw', '.bed.gz'),
        tbi=data_path + reference_path + gerp_url.split('/')[-1].replace('.bw', '.bed.gz.tbi')
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip {input}; tabix {output.bed}"

# chrom start end lowerLimit upperLimit average
rule download_phastcons:
    output:
        bed=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz', '.bed.gz'),
        tbi=data_path + reference_path + phastCons_url.split('/')[-1].replace('.txt.gz', '.bed.gz.tbi')
    params:
        url = phastCons_url,
        data_path = data_path,
        reference_path = reference_path,
        tmp_file= data_path + reference_path + phastCons_url.split('/')[-1]
    conda:
        "../envs/tabix.yaml"
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}; zcat {params.tmp_file} | "
        "awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $2, $3, $4, $10, $10 + $11, $13 / $12}}' | "
        "grep -E -w '^chr[0-9][0-9]?' | sort -k1,1 -k2,2n | bgzip > {output.bed}; "
        "tabix {output.bed}; rm {params.tmp_file}"

# chrom start end lowerLimit upperLimit average
use rule download_phastcons as download_phylop with:
    output:
        bed=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz', '.bed.gz'),
        tbi=data_path + reference_path + phyloP_url.split('/')[-1].replace('.txt.gz', '.bed.gz.tbi')
    params:
        url = phyloP_url,
        data_path = data_path,
        reference_path = reference_path,
        tmp_file= data_path + reference_path + phyloP_url.split('/')[-1]

rule download_recomb_rate:
    output:
        temp(data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed'))
    params:
        url = recomb_rate_url,
        data_path = data_path,
        reference_path = reference_path,
        prefix = recomb_rate_url.split('/')[-1].split(".bw")[0]
    conda:
        "../envs/bigwigtowig.yaml"
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}; "
        "bigWigToWig {params.data_path}{params.reference_path}{params.prefix}.bw "
        "{params.data_path}{params.reference_path}{params.prefix}.wig; "
        "wig2bed < {params.data_path}{params.reference_path}{params.prefix}.wig | "
        "grep -E -w '^chr[0-9][0-9]?' | sort -k1,1 -k2,2n > {output}; rm {params.data_path}{params.reference_path}{params.prefix}.bw; "
        "rm {params.data_path}{params.reference_path}{params.prefix}.wig"

use rule compress_and_index_gerp_scores as compress_and_index_rec_rate with:
    input:
        data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw','.bed')
    output:
        bed=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz'),
        tbi=data_path + reference_path + recomb_rate_url.split('/')[-1].replace('.bw', '.bed.gz.tbi')
#
rule download_bstat:
    output:
        temp(data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0] + '.bed'),
    params:
        url = bstat_url
    shell:
        "wget -qO- {params.url} | grep -E -w '^chr[0-9]+' | sort -k1,1 -k2,2n > {output}"

use rule lift_1kgp_mask_from_hg19_to_hg38 as lift_bstats_from_hg19_to_hg38 with:
    input:
        bed=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0] + '.bed',
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz"
    output:
        temp(data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed')

rule compress_and_index_bstat:
    input:
        data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed'
    output:
        bed=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed.gz',
        tbi=data_path+ reference_path + bstat_url.split('/')[-1].split('.txt')[0].replace('hg19', 'hg38') + '.bed.gz.tbi'
    conda:
        "../envs/tabix.yaml"
    shell:
        'sort -k1,1 -k2,2n {input} | bgzip > {output.bed}; tabix {output.bed}'

rule download_gnomad:
    output:
        temp(data_path + reference_path + gnomad_url.split('/')[-1]),
        temp(data_path + reference_path + gnomad_url.split('/')[-1] + '.tbi')
    params:
        url = gnomad_url,
        data_path = data_path,
        reference_path = reference_path
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}; "
        "wget -q -P {params.data_path}{params.reference_path} {params.url}.tbi"

use rule download_gnomad as download_clinvar with:
    output:
        temp(data_path + reference_path + clinvar_url.split('/')[-1]),
        temp(data_path + reference_path + clinvar_url.split('/')[-1] + '.tbi')
    params:
        url = clinvar_url,
        data_path= data_path,
        reference_path= reference_path

rule download_gtex:
    output:
        temp(data_path + reference_path + gtex_url.split('/')[-1])
    params:
        url = gtex_url,
        data_path= data_path,
        reference_path= reference_path
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}"

checkpoint unpack_gtex:
    input:
        data_path + reference_path + gtex_url.split('/')[-1]
    output:
        directory(data_path + reference_path + gtex_url.split('/')[-1].replace('.tar', ''))
    params:
        output_dir = data_path + reference_path
    shell:
        "tar xf {input} -C {params.output_dir}"


def get_number_of_tissues_gtex(wildcards):
    checkpoints.unpack_gtex.get(**wildcards)
    tissues = glob_wildcards(data_path + reference_path +  gtex_url.split('/')[-1].replace('.tar', '') +
                             "/{tissue}.v8.independent_eqtls.txt.gz").tissue
    return len(tissues)

# do bonferroni correction --> pval_nominal_threshold/n where n is the number of tissue
# chrom start end pval_nominal pval_thresh qval gene_id tissue
rule get_significant_eqtls:
    input:
        variants=data_path + reference_path + gtex_url.split('/')[-1].replace('.tar','') +
                 "/{tissue}.v8.independent_eqtls.txt.gz"
    output:
        temp(data_path + reference_path + "{tissue}.v8.independent_genome_wide_significant_eqtls.txt.gz")
    params:
        tissue = '{tissue}',
        n_tissues = get_number_of_tissues_gtex
    shell:
        "zcat {input.variants} | sed 's/_/\t/g' | awk -v tissue={params.tissue} -v n_tissues={params.n_tissues} "
        "-F '\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($17 < 5e-8 / n_tissues) print $7, $8 - 1, $8,  $1, $17, tissue}}' > {output}"

def get_significant_eqtls_all_tissues(wildcards):
    checkpoints.unpack_gtex.get(**wildcards)
    tissues = glob_wildcards(data_path + reference_path +  gtex_url.split('/')[-1].replace('.tar', '') +
                             "/{tissue}.v8.independent_eqtls.txt.gz").tissue
    return expand(data_path + reference_path +  "{tissue}.v8.independent_genome_wide_significant_eqtls.txt.gz", tissue=tissues)

rule aggregate_significant_eqtls:
    input:
        get_significant_eqtls_all_tissues
    output:
        data_path + reference_path + "GTEx_significant_independent_eqtls_all_tissues.bed"
    resources:
        mem_mb=32 * 1000
    shell:
        "cat {input} | sort -k1,1 -k2,2n > {output}"

# chrom start end type score
rule download_encode_annotation:
    output:
        bed=temp(data_path + reference_path + 'tmp_' + encode_annotation_url.split('/')[-1].replace('.txt.gz','.bed'))
    params:
        url = encode_annotation_url,
        data_path = data_path,
        reference_path = reference_path,
        tmp_file= data_path + reference_path + encode_annotation_url.split('/')[-1]
    shell:
        "wget -q -P {params.data_path}{params.reference_path} {params.url}; zcat {params.tmp_file} | "
        "awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{if ($5 != \"R\" && $5 != \"T\") print $2, $3, $4, $5, $6}}' | "
        "grep -E -w '^chr[0-9][0-9]?' | sort -k1,1 -k2,2n > {output.bed}; rm {params.tmp_file}"


use rule lift_1kgp_mask_from_hg19_to_hg38 as lift_encode_annotation_from_hg19_to_hg38 with:
    input:
        bed=data_path + reference_path + 'tmp_' + encode_annotation_url.split('/')[-1].replace('.txt.gz','.bed'),
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz"
    output:
        bed=temp(data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','.bed'))

use rule compress_and_index_bstat as compress_and_index_encode_anno with:
    input:
        data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','.bed')
    output:
        bed=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','.bed.gz'),
        tbi=data_path + reference_path + encode_annotation_url.split('/')[-1].replace('.txt.gz','.bed.gz.tbi')

# download chromosome sizes
rule download_chromosome_sizes:
    output:
        data_path + reference_path + hg38_chrom_sizes_url.split('/')[-1] + '.bed'
    params:
        url = hg38_chrom_sizes_url
    shell:
        "wget -qO- {params.url} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1, 0, $2}}' | grep -w 'chr[0-9]*' | "
        "sort -k1,1 -k2,2n >{output}"

rule make_genomefile:
    input:
        data_path + reference_path + hg38_chrom_sizes_url.split('/')[-1] + '.bed'
    output:
        data_path + reference_path + 'genomefile_hg38.bed'
    shell:
        "cut -f1,3 {input} > {output}"

rule window_genome:
    input:
        data_path + reference_path + 'genomefile_hg38.bed'
    output:
        data_path + reference_path + "hg38_windowed_w_{windowsize}_s_{stepsize}.bed"
    params:
        bedtools = bedtools_path
    shell:
        "{params.bedtools} makewindows -g {input} -w {wildcards.windowsize} -s {wildcards.stepsize} > {output}"

# human reference genome
rule download_human_reference_genome:
    output:
        data_path + reference_path + hg38_fasta_url.split('/')[-1].replace('.gz', '')
    params:
        url = hg38_fasta_url,
        path = data_path + reference_path
    shell:
        "wget -q -P {params.path} {params.url}; gunzip {output}.gz"

rule download_cytobands:
    output:
        data_path + reference_path + 'hg38_cytobands.txt'
    params:
        url = hg38_cytobands_url
    shell:
        'wget -q -O {output}.gz {params.url}; gunzip {output}.gz'

rule get_chromosome_arms:
    input:
        data_path + reference_path + 'hg38_cytobands.txt'
    output:
        data_path + reference_path + 'hg38_chromosome_arms.bed'
    run:
        cytobands = pd.read_csv(input[0],sep='\t',names=["chrom", "start", "end", "name", "gieStain"])
        cytobands.chrom = [chrom.replace('chr','') for chrom in cytobands.chrom]
        cytobands = cytobands[np.isin(cytobands.chrom.values,np.arange(1, 23).astype(str))]
        cytobands.name = [name[0] for name in cytobands.name]
        cytobands.chrom = cytobands.chrom.astype(int)
        cytobands.start = cytobands.start.astype(int)
        cytobands.end = cytobands.end.astype(int)
        chromosome_arms = cytobands.groupby(["chrom", 'name']).agg({'start': np.min, 'end': np.max}).reset_index()
        chromosome_arms = chromosome_arms.loc[:, ['chrom', 'start', 'end', 'name']].rename({'name': 'arm'},axis=1)
        chromosome_arms.sort_values(['chrom', 'start', 'end'],inplace=True)
        chromosome_arms.to_csv(output[0], sep='\t', index=False, header=False)

def get_url_whole_genome_axt(wildcards):
    # return whole_genome_alignments_url.format(species_upper=wildcards.species[0].upper() + wildcards.species[1:],
    #     species=wildcards.species, chr=wildcards.chr)
    return whole_genome_alignments_url.format(species_upper=wildcards.species[0].upper() + wildcards.species[1:],
        species=wildcards.species)

# download axt file of human-primates whole genome alignments
rule download_whole_genome_alignments_axt:
    output:
        temp(data_path + whole_genome_alignments_path + whole_genome_alignments_url.split('/')[-1])
    params:
        data_path = data_path,
        whole_genome_alignments_path = whole_genome_alignments_path,
        url = get_url_whole_genome_axt
    shell:
        "wget -q -P {params.data_path}{params.whole_genome_alignments_path} {params.url}"

# download Li Heng's seqbility
rule download_and_install_seqbility:
    output:
        directory("seqbility")
    params:
        url = seqbility_url
    shell:
        "wget -q {params.url}; tar -xf {output}-20091110.tar.bz2; mv {output}-20091110 {output}; cd {output}; make"

# extract sample IDs for populations of interest
rule extract_sample_ids_of_population:
    input:
       samples = data_path + "1000G_phase3/integrated_call_samples_v3.20130502.ALL.panel",
    output:
       data_path + "{population}_sample_ids.txt"
    wildcard_constraints:
        population='|'.join([pop for pop in populations if pop != 'AA' and pop != 'AOUAFR'
                             and pop != 'AOUEUR' and pop != 'AOUNA'])
    shell:
        "if [ '{wildcards.population}' != 'ALL' ]; then "
        "awk -F '\\t' '{{if ($2 == \"{wildcards.population}\") print $1}}' {input} > {output}; "
        "else "
        "tail -n+2 {input} | cut -f1 > {output}; "
        "fi"

def get_pops_of_superpop(wildcards):
    return [f"{data_path}{pop}_sample_ids.txt" for pop, superpop in pop_superpop_mapping.items()
            if superpop == wildcards.superpopulation and pop != 'AA' and pop != 'AOUAFR'
            and pop != 'AOUEUR' and pop != 'AOUNA']

rule extract_sample_ids_of_superpopulation:
    input:
        get_pops_of_superpop
    output:
       data_path + "{superpopulation}_sample_ids.txt"
    wildcard_constraints:
        superpopulation='AFR|EUR|EAS'
    shell:
        "cat {input} > {output}"

def get_input_create_reference_panel(wildcards):
    afr = ["ESN", "GWD", "LWK", "MSL", "YRI"]
    eur = ["CEU", "FIN", "GBR", "IBS", "TSI"]
    eas = ["CDX", "CHS", "CHB", "KHV", "JPT"]
    return [data_path + f"{pop}_sample_ids.txt" for populations in [afr, eur, eas] for pop in populations]

rule create_reference_panel:
    input:
        get_input_create_reference_panel
    output:
        data_path + "1000G_phase3/reference_panel.txt"
    params:
        afr = "ESN GWD LWK MSL YRI",
        eur = "CEU FIN GBR IBS TSI",
        eas = "CDX CHS CHB KHV JPT",
        data_path = data_path
    shell:
        "for pop in {params.afr}; do awk 'BEGIN{{OFS=\"\\t\"}}{{print $1, \"AFR\"}}' {params.data_path}${{pop}}_sample_ids.txt >> {output}; done; "
        "for pop in {params.eur}; do awk 'BEGIN{{OFS=\"\\t\"}}{{print $1, \"EUR\"}}' {params.data_path}${{pop}}_sample_ids.txt >> {output}; done; "
        "for pop in {params.eas}; do awk 'BEGIN{{OFS=\"\\t\"}}{{print $1, \"EAS\"}}' {params.data_path}${{pop}}_sample_ids.txt >> {output}; done; "

# created based on HapMap data
rule download_genetic_map_plink:
    output:
        expand(data_path + genetic_map_path + "plink.chr{chr}.GRCh38.map", chr=chromosomes)
    params:
        output_path = data_path + genetic_map_path,
        url = genetic_map_plink_url,
        zipped_file = genetic_map_plink_url.split('/')[-1]
    shell:
        "mkdir -p {params.output_path}; wget -q {params.url}; unzip -d {params.output_path} {params.zipped_file}; "
        "rm {params.zipped_file}"

rule format_genetic_map:
    input:
        data_path + genetic_map_path + "plink.chr{chr}.GRCh38.map"
    output:
        data_path + genetic_map_path + "plink.formatted.chr{chr}.GRCh38.map"
    shell:
        "sed -e 's/^/chr/' {input} > {output}"

rule download_genetic_map:
    output:
        expand(data_path + genetic_map_path + "genetic_map_Hg38_chr{chr}.txt",chr=chromosomes)
    params:
        output_path = data_path + genetic_map_path,
        url = genetic_map_url,
        compressed_file = data_path + genetic_map_path + genetic_map_url.split('/')[-1],
        filenames="genetic_map_Hg38_chr*.txt"
    shell:
        "mkdir -p {params.output_path}; wget -q -P {params.output_path} {params.url}; tar xzvf {params.compressed_file};"
        "rm {params.compressed_file}; mv {params.filenames} {params.output_path}"

rule get_individuals_from_reference_populations:
    input:
        pops = expand(data_path + "{population}_sample_ids.txt",population=[pop for pop in populations if pop != 'AA'
                                                                            and pop != 'AOUAFR'
                                                                            and pop != 'AOUEUR' and pop != 'AOUNA'])
    output:
        reference_pops = data_path + "1000G_phase3/individuals_in_reference_populations_1kgp.txt"
    shell:
        "cat {input.pops} > {output.reference_pops}; "

# extract reference populations from 1KGP
rule extract_ref_pops_from_1KGP:
    input:
        vcf = data_path + "1000G_phase3/ALL.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        pops = data_path + "1000G_phase3/individuals_in_reference_populations_1kgp.txt"
    output:
        vcf=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    params:
        bcftools = bcftools_path
    threads: 8
    shell:
        "{params.bcftools} view -S {input.pops} -M 2 -v snps {input.vcf} | {params.bcftools} norm -d exact -Oz -o {output.vcf} --threads {threads}"

rule index_ref_pops_1kgp:
    input:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
    output:
        data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
    conda:
        '../envs/tabix.yaml'
    shell:
        "tabix {input}"

# convert 1KGP VCF to plink file format
rule vcf_to_plink_1kgp:
    input:
        vcf=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        fasta=data_path + reference_path + "hg38.fa"
    output:
        temp(multiext(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes",
                      ".bed", ".bim", ".fam"))
    params:
        base = data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes"
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    threads: 8
    shell:
        "plink2 --vcf {input.vcf} --snps-only --chr 1-22 --allow-extra-chr --max-alleles 2 --memory {resources.mem_mb} "
        "--set-all-var-ids @:#\$1:\$2 --make-bed --out {params.base} --threads {threads} --ref-from-fa {input.fasta}"

rule update_1kgp_fam:
    input:
        plink= multiext(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes",
                        ".bed", ".bim", ".fam"),
        panel=data_path + "1000G_phase3/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        multiext(data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.updated_ids",
                      ".bed", ".bim", ".fam")
    params:
        input_base=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes",
        output_base=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.updated_ids"
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    threads: 8
    shell:
        "cut -f2 {params.input_base}.fam | xargs -I {{}} grep {{}} {input.panel} | awk -F '\\t' '{{print $2,$1}}' | "
        "paste <(cut -f1,2 {params.input_base}.fam) - > {params.output_base}.fam.tmp; "
        "plink2 --bfile {params.input_base} --update-ids {params.output_base}.fam.tmp --make-bed "
        "--out {params.output_base} --threads {threads} --memory {resources.mem_mb}; "
        "rm {params.output_base}.fam.tmp"


rule get_hgdp_population_sample_ids:
     output:
        maya=data_path + "Maya_sample_ids.txt",
        pima=data_path + "Pima_sample_ids.txt"
     params:
         maya=hgdp_mayan_samples,
         pima=hgdp_pima_samples
     shell:
         "for sample in {params.maya}; do echo ${{sample}} >> {output.maya}; done; "
         "for sample in {params.pima}; do echo ${{sample}} >> {output.pima}; done; "

#
rule vcf_to_plink_hgdp:
    input:
        vcf=data_path + "hgdp/hgdp_wgs.20190516.full.chr{chr}.vcf.gz",
        maya=data_path + "Maya_sample_ids.txt",
        pima=data_path + "Pima_sample_ids.txt"
    output:
        temp(multiext(data_path + "hgdp/hgdp_wgs.20190516.na.chr{chr}", ".bed",".bim",".fam"))
    params:
        base=data_path + "hgdp/hgdp_wgs.20190516.na.chr{chr}"
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    threads: 8
    shell:
        "plink2 --vcf {input.vcf} --snps-only --allow-extra-chr --max-alleles 2 --memory {resources.mem_mb} "
        "--set-all-var-ids @:#\$1:\$2 --make-bed --out {params.base} --threads {threads} "
        "--keep <(cat {input.maya} {input.pima})"

rule update_ids_hgdp:
    input:
        multiext(data_path + "hgdp/hgdp_wgs.20190516.na.chr{chr}",".bed",".fam",".bim")
    output:
        multiext(data_path + "hgdp/hgdp_wgs.20190516.native_americans.chr{chr}", ".bed", ".fam", ".bim")
    conda:
        "../envs/plink.yaml"
    threads: 16
    params:
        input_prefix=data_path + "hgdp/hgdp_wgs.20190516.na.chr{chr}",
        output_prefix=data_path + "hgdp/hgdp_wgs.20190516.native_americans.chr{chr}",
        hgdp_mayan_samples=hgdp_mayan_samples,
        hgdp_pima_samples=hgdp_pima_samples
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "for sample in {params.hgdp_mayan_samples}; do echo -e \"0\\t${{sample}}\\tMayan\\t${{sample}}\"; done "
        "> native_american_samples_updated_ids_{wildcards.chr}.txt; "
        "for sample in {params.hgdp_pima_samples}; do echo -e \"0\\t${{sample}}\\tPima\\t${{sample}}\"; done "
        ">> native_american_samples_updated_ids_{wildcards.chr}.txt; "
        "plink2 --bfile {params.input_prefix} --update-ids native_american_samples_updated_ids_{wildcards.chr}.txt "
        "--memory {resources.mem_mb} --make-bed "
        "--out {params.output_prefix} --threads {threads}; "
        "rm native_american_samples_updated_ids_{wildcards.chr}.txt;"

#
# rule index_target_data_vcf:
#     input:
#         target_data_vcf
#     output:
#         target_data_vcf.replace('.vcf.gz', '.vcf.gz.tbi')
#     conda:
#         "../envs/tabix.yaml"
#     shell:
#         "tabix {input}"

# update variant IDs in target plink files so that they match
# rule vcf_to_plink_target_data:
#     input:
#         vcf = target_data_vcf,
#         tbi = target_data_vcf.replace('.vcf.gz','.vcf.gz.tbi'),
#         fasta=data_path + reference_path + "hg38.fa"
#     output:
#         temp(multiext(target_data_vcf.replace('.vcf.gz', '_custom_variant_ids'), ".bim", ".bed", ".fam"))
#     params:
#         output_base = target_data_vcf.replace('.vcf.gz', '_custom_variant_ids')
#     threads: 8
#     conda:
#         "../envs/plink.yaml"
#     resources:
#         mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
#     shell:
#         "plink2 --vcf {input.vcf} --snps-only --allow-extra-chr --max-alleles 2 "
#         "--set-all-var-ids @:#\$1:\$2 --make-bed --out {params.output_base} --threads {threads} "
#         "--memory {resources.mem_mb} --ref-from-fa {input.fasta}"

#
# rule get_modern_vcf:
#     input:
#         plink=multiext(data_path + merged_datasets_path + "1KG_HGDP_chr{chr}", '.bed', '.bim', '.fam'),
#         fasta=data_path + reference_path + "hg38.fa"
#     output:
#         data_path + merged_datasets_path + "1KG_HGDP_chr{chr}.vcf.gz"
#     params:
#         base=data_path + merged_datasets_path + "1KG_HGDP_chr{chr}"
#     conda:
#         "../envs/plink.yaml"
#     threads: 8
#     shell:
#         "plink2 --bfile {params.base} --export vcf bgz id-paste=iid --out {params.base} --threads {threads} "
#         "--ref-from-fa {input.fasta}"

# # merge with Native American samples from HGDP
# rule get_shared_variants_target:
#     input:
#         target = target_data_vcf.replace('.vcf.gz', '_custom_variant_ids.bim'),
#         kgp_hgdp = data_path + merged_datasets_path + "1KG_HGDP_chr{chr}.bim"
#     output:
#         temp("shared_variants_chr{chr}_target_kgp_hgdp.txt")
#     wildcard_constraints:
#         chr="[0-9]+"
#     shell:
#         "comm -12 <(sort <(cut -f2 {input.kgp_hgdp})) <(sort <(cut -f2 {input.target})) > {output}"

# use rule extract_shared_variants_hgdp_data as extract_shared_variants_kgp_hgdp_data with:
#     input:
#         shared = "shared_variants_chr{chr}_target_kgp_hgdp.txt",
#         bed = multiext(data_path + merged_datasets_path + "1KG_HGDP_chr{chr}", ".bed", ".bim", ".fam")
#     output:
#         temp(multiext(data_path + merged_datasets_path + "1KG_HGDP_chr{chr}_shared_variants",
#             ".bed", ".bim", ".fam"))
#     params:
#         input_base = data_path + merged_datasets_path + "1KG_HGDP_chr{chr}",
#         output_base = data_path + merged_datasets_path + "1KG_HGDP_chr{chr}_shared_variants"

# use rule extract_shared_variants_hgdp_data as extract_shared_variants_target_data with:
#     input:
#         shared = "shared_variants_chr{chr}_target_kgp_hgdp.txt",
#         bed = multiext(target_data_vcf.replace('.vcf.gz', '_custom_variant_ids'), ".bed", ".bim", ".fam")
#     output:
#         temp(multiext(target_data_vcf.replace('.vcf.gz', '_custom_variant_ids_shared_variants'), ".bed", ".bim", ".fam"))
#     params:
#         input_base = target_data_vcf.replace('.vcf.gz', '_custom_variant_ids'),
#         output_base = target_data_vcf.replace('.vcf.gz', '_custom_variant_ids_shared_variants')

# use rule merge_datasets as merge_datasets_with_target with:
#     input:
#         reference = multiext(data_path + merged_datasets_path + "1KG_HGDP_chr{chr}_shared_variants",
#                              ".bed", ".bim", ".fam"),
#         target = multiext(target_data_vcf.replace('.vcf.gz', '_custom_variant_ids_shared_variants'),
#                           ".bed", ".bim", ".fam")
#     output:
#         temp(multiext(data_path + merged_datasets_path + "1KG_HGDP_target_chr{chr}", ".bed", ".bim", ".fam"))
#     params:
#         output_base = data_path + merged_datasets_path + "1KG_HGDP_target_chr{chr}",
#         reference_base = data_path + merged_datasets_path + "1KG_HGDP_chr{chr}_shared_variants",
#         target_base = target_data_vcf.replace('.vcf.gz', '_custom_variant_ids_shared_variants'),
#         maf = maf,
#         geno = geno
#     conda:
#         "../envs/plink.yaml"
#     resources:
#         mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000

rule merged_dataset_to_plink:
    input:
        data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.vcf.gz'
    output:
        temp(multiext(data_path  + merged_datasets_path + 'tmp_kgp_target_individuals_chr{chr}', ".bed", ".bim", ".fam"))
    conda:
        "../envs/plink.yaml"
    params:
        prefix=data_path  + merged_datasets_path + 'tmp_kgp_target_individuals_chr{chr}',
        maf = maf,
        geno = geno
    threads: 1
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --vcf {input} --make-bed --out {params.prefix} --set-all-var-ids @:#\$1:\$2 --threads {threads} "
        "--maf {params.maf} --geno {params.geno} --memory {resources.mem_mb}"

rule update_ids_merged_dataset_plink:
    input:
        plink=multiext(data_path + merged_datasets_path  + 'tmp_kgp_target_individuals_chr{chr}',".bed",".bim",".fam"),
        ref_panel=data_path + "1000G_phase3/reference_panel.txt",
        indv_pop =data_path + "1000G_phase3/integrated_call_samples_v3.20130502.ALL.panel"
    output:
        plink=multiext(data_path + merged_datasets_path  + 'kgp_target_individuals_chr{chr}',".bed",".bim",".fam"),
        updated_ids=temp("updated_ids_{chr}.txt"),
        updated_ids_tmp=temp("tmp_updated_ids_{chr}.txt")
    params:
        input_prefix=data_path + merged_datasets_path  + 'tmp_kgp_target_individuals_chr{chr}',
        output_prefix=data_path + merged_datasets_path + 'kgp_target_individuals_chr{chr}'
    conda:
        "../envs/plink.yaml"
    threads: 1
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "paste <(cut -f1 {input.ref_panel} | grep -f - {params.input_prefix}.fam | cut -f1,2 | sort -k2) "
        "<(cut -f1 {input.ref_panel} | grep -f - {input.indv_pop} | awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $2, $1}}' | sort -k2) "
        "> {output.updated_ids_tmp}; "
        "cut -f1 {input.ref_panel} | grep -v -f - {params.input_prefix}.fam | "
        "awk -F '\\t' 'BEGIN{{OFS=\"\\t\"}}{{print $1, $2, $1, $2}}' >> {output.updated_ids_tmp}; "
        "cut -f2  {params.input_prefix}.fam | xargs -I {{}} grep {{}} {output.updated_ids_tmp} > {output.updated_ids}; "
        "plink2 --bfile {params.input_prefix} --update-ids {output.updated_ids} --make-bed --out {params.output_prefix} "
        "--threads {threads} --memory {resources.mem_mb}"


# merge HGDP data with other datasets
rule get_shared_variants:
    input:
        hgdp=data_path + "hgdp/hgdp_wgs.20190516.native_americans.chr{chr}.bim",
        kgp = data_path  + merged_datasets_path + 'kgp_target_individuals_chr{chr}.bim'
    output:
        temp("shared_variants_chr{chr}.txt")
    wildcard_constraints:
        chr="[0-9]*"
    shell:
        "comm -12 <(sort <(cut -f2 {input.hgdp})) <(sort <(cut -f2 {input.kgp})) > {output}"


# # extract the shared variants
rule extract_shared_variants_hgdp_data:
    input:
        shared = 'shared_variants_chr{chr}.txt',
        plink=multiext(data_path + "hgdp/hgdp_wgs.20190516.native_americans.chr{chr}", ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + "hgdp/hgdp_wgs.20190516.native_americans_shared_variants.chr{chr}",".bed",".bim",".fam"))
    params:
        input_base = data_path + "hgdp/hgdp_wgs.20190516.native_americans.chr{chr}",
        output_base= data_path + "hgdp/hgdp_wgs.20190516.native_americans_shared_variants.chr{chr}"
    threads: 8
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    wildcard_constraints:
        chr="[0-9]+"
    shell:
        "plink2 --bfile {params.input_base} --extract {input.shared} --make-bed --out {params.output_base} "
        "--threads {threads} --memory {resources.mem_mb}"
#
use rule extract_shared_variants_hgdp_data as extract_shared_variants_kgp_target with:
    input:
        shared='shared_variants_chr{chr}.txt',
        plink=multiext(data_path  + merged_datasets_path + 'kgp_target_individuals_chr{chr}',
                       ".bed", ".bim", ".fam")
    output:
        temp(multiext(data_path + merged_datasets_path + 'hgdp_kgp_target_individuals_chr{chr}',
            ".bed", ".bim", ".fam"))
    wildcard_constraints:
        chr="[0-9]+"
    params:
        input_base = data_path + merged_datasets_path  + 'kgp_target_individuals_chr{chr}',
        output_base = data_path + merged_datasets_path  + 'hgdp_kgp_target_individuals_chr{chr}'
#
# merge the two datasets using only shared variants
rule merge_datasets:
    input:
        reference = multiext(data_path + merged_datasets_path  + 'kgp_target_individuals_chr{chr}',
                             ".bed", ".bim", ".fam"),
        target=multiext(data_path + "hgdp/hgdp_wgs.20190516.native_americans_shared_variants.chr{chr}",
                        ".bed", ".bim", ".fam"),
    output:
        multiext(data_path + merged_datasets_path + 'reference_target_individuals_chr{chr}', ".bed", ".bim", ".fam")
    params:
        output_base = data_path + merged_datasets_path + 'reference_target_individuals_chr{chr}',
        reference_base = data_path + merged_datasets_path + 'kgp_target_individuals_chr{chr}',
        target_base= data_path + "hgdp/hgdp_wgs.20190516.native_americans_shared_variants.chr{chr}",
        maf = maf,
        geno = geno
    conda:
        "../envs/plink.yaml"
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    threads: 1
    retries: 2
    shell:
        "if [[ -f {params.output_base}_trial1-merge.missnp ]]; "  # second trial --> flip SNPs once
        "then "
            "plink --bfile {params.target_base} --flip {params.output_base}_trial1-merge.missnp --make-bed --out "
            "{params.target_base}_tmp --threads {threads}; "  # flip once
            "rm {params.output_base}_trial1-merge.missnp; "
            "plink --bfile {params.reference_base} --bmerge {params.target_base}_tmp --make-bed --geno {params.geno} "
            "--maf {params.maf} --out {params.output_base}_trial2 --threads {threads} --memory {resources.mem_mb}; "
            "mv {params.output_base}_trial2.bed {params.output_base}.bed; "
            "mv {params.output_base}_trial2.bim {params.output_base}.bim; "
            "mv {params.output_base}_trial2.fam {params.output_base}.fam; "
        "elif [[ -f {params.output_base}_trial2-merge.missnp ]]; "  # last trial --> excluding SNPS that caused problems
        "then "
            "plink --bfile {params.reference_base} --exclude {params.output_base}_trial2-merge.missnp --make-bed "
            "--out {params.reference_base}_tmp --threads {threads}; "
            "plink --bfile {params.target_base}_tmp --exclude {params.output_base}_trial2-merge.missnp --make-bed "
            "--out {params.target_base}_tmp1 --threads {threads}; "
            "plink --bfile {params.reference_base}_tmp --bmerge {params.target_base}_tmp1 --make-bed --geno {params.geno} "
            "--maf {params.maf} --out {params.output_base} --threads {threads} --memory {resources.mem_mb}; "
            "rm {params.reference_base}_tmp.*; "
            "rm {params.target_base}_tmp.*; "
            "rm {params.target_base}_tmp1.*; "
            "rm {params.output_base}_trial2-merge.missnp; "
        "else "  # first pass
            "plink --bfile {params.reference_base} --bmerge {params.target_base} --make-bed --geno {params.geno} "
            "--maf {params.maf} --out {params.output_base}_trial1 --threads {threads} --memory {resources.mem_mb}; "
            "mv {params.output_base}_trial1.bed {params.output_base}.bed; "
            "mv {params.output_base}_trial1.bim {params.output_base}.bim; "
            "mv {params.output_base}_trial1.fam {params.output_base}.fam; "
        "fi"


rule merge_chromosomes_1kgp_hgdp_target:
    input:
        bim=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bim', chr=chromosomes),
        fam=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.fam', chr=chromosomes),
        bed=expand(data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.bed', chr=chromosomes),
        kin=related_samples
    output:
        temp(multiext(data_path + merged_datasets_path  + 'reference_target_individuals_ALL', ".bed", ".bim", ".fam")),
        bedfiles=temp("bed_files_1kgp_hgdp_target.txt")
    params:
        output_base= data_path + merged_datasets_path  + 'reference_target_individuals_ALL'
    threads: 8
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    conda:
        "../envs/plink.yaml"
    shell:
        "for bed in {input.bed}; do echo $bed | sed 's/.bed//' >> {output.bedfiles}; done; "
        "plink --merge-list {output.bedfiles} --make-bed --allow-no-sex --out {params.output_base} --threads {threads} "
        "--memory {resources.mem_mb}"

rule filter_kin_reference_panel:
    input:
        kin = related_samples,
        ref = data_path + "1000G_phase3/reference_panel.txt"
    output:
        data_path + "1000G_phase3/reference_panel_filtered.txt"
    shell:
        "grep -v -f {input.kin} {input.ref} > {output}"

rule merge_reference_and_target_vcf:
    input:
        ref=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz",
        ref_tbi=data_path + "1000G_phase3/REF.chr{chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi",
        target=data_path + 'target_individuals_chr{chr}.vcf.gz',
        target_tbi=data_path + 'target_individuals_chr{chr}.vcf.gz.tbi',
    output:
        data_path + merged_datasets_path  + 'reference_target_individuals_chr{chr}.vcf.gz'
    params:
        bcftools = bcftools_path
    shell:
        "{params.bcftools} merge {input.ref} {input.target} | "
        "{params.bcftools} view -v snps -M 2 -i 'F_MISSING==0.0' | "
        "{params.bcftools} norm -d exact -Oz -o {output}"

# # LD prune merged dataset
rule ld_prune_combined_dataset_1kgp_hgdp_target:
    input:
        multiext(data_path + merged_datasets_path  + 'reference_target_individuals_ALL', ".bed", ".bim", ".fam")
    output:
        multiext(data_path + merged_datasets_path  + 'reference_target_individuals_ALL', ".prune.in", ".prune.out")
    params:
        output_base = data_path + merged_datasets_path  + 'reference_target_individuals_ALL',
        input_base = data_path + merged_datasets_path + 'reference_target_individuals_ALL'
    conda:
        "../envs/plink.yaml"
    threads: 16
    resources:
        mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --bfile {params.input_base} --indep-pairwise 200 kb 1 0.1 --out {params.output_base} --threads {threads} "
        "--memory {resources.mem_mb}"

# filter out variants in LD
rule filter_ld_pruned_variants_1kgp_hgdp:
    input:
        multiext(data_path + merged_datasets_path  + 'reference_target_individuals_ALL',".bim",".bed",".fam"),
        prunein=data_path + merged_datasets_path  + 'reference_target_individuals_ALL.prune.in'
    output:
        multiext(data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned', ".bed", ".bim", ".fam")
    params:
        input_base = data_path + merged_datasets_path + 'reference_target_individuals_ALL',
        output_base = data_path + merged_datasets_path  + 'reference_target_individuals_ALL_ld_pruned'
    conda:
        "../envs/plink.yaml"
    threads: 16
    resources:
        mem_mb=lambda wildcards, threads: threads * mem_gb_cpu * 1000
    shell:
        "plink2 --bfile {params.input_base} --extract {input.prunein} --make-bed --out {params.output_base}  --threads {threads} "
        "--memory {resources.mem_mb}"

# rule extract_shared_variants_target:
#     input:
#         plink=expand(data_path + 'target_individuals_chr{chr}{ext}', chr=chromosomes, ext=['.bed', '.bim', '.fam']),
#         shared=data_path + merged_datasets_path + "shared_variants_reference_and_target.txt"
#     output:
#         temp(multiext(data_path + merged_datasets_path + "target_dataset_shared_variants", ".bed", ".bim", ".fam")),
#         bedfiles=temp("bed_files_target.txt")
#     params:
#         input_prefix = expand(data_path + 'target_individuals_chr{chr}', chr=chromosomes),
#         output_prefix = data_path + merged_datasets_path + "target_dataset_shared_variants"
#     conda:
#         "../envs/plink.yaml"
#     threads: 16
#     resources:
#         mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000
#     shell:
#         "for plink in {params.input_prefix}; do echo ${{plink}} >> {output.bedfiles}; done; "
#         "plink --merge-list {output.bedfiles} --extract {input.shared} --make-bed --out {params.output_prefix} "
#         "--threads {threads} --memory {resources.mem_mb}"

# use rule merge_datasets as merge_datasets_with_target with:
#     input:
#         reference = (multiext(data_path + merged_datasets_path + "1KG_HGDP_ALL_shared_variants_w_target", ".bim", ".bed", ".fam")),
#         target = multiext(data_path + merged_datasets_path + "target_dataset_shared_variants", ".bed", ".bim", ".fam")
#     output:
#         multiext(data_path + merged_datasets_path + "1KG_HGDP_target_ALL", ".bed", ".bim", ".fam")
#     params:
#         output_base = data_path + merged_datasets_path + "1KG_HGDP_target_ALL",
#         reference_base = data_path + merged_datasets_path + "1KG_HGDP_ALL_shared_variants_w_target",
#         target_base = data_path + merged_datasets_path + "target_dataset_shared_variants",
#         maf = maf,
#         geno = geno
#     conda:
#         "../envs/plink.yaml"
#     resources:
#         mem_mb = lambda wildcards, threads: threads * mem_gb_cpu * 1000

rule download_boosting_scores_complete:
    output:
        eur=temp(data_path + reference_path + "tmp_EUR_complete_boosting_scores.bed.gz"),
        eur_recent=temp(data_path+ reference_path + "tmp_EUR_complete_recent_boosting_scores.bed.gz"),
        eur_ancient=temp(data_path+ reference_path + "tmp_EUR_complete_ancient_boosting_scores.bed.gz"),
        afr=temp(data_path + reference_path + "tmp_AFR_complete_boosting_scores.bed.gz"),
        afr_recent=temp(data_path + reference_path + "tmp_AFR_complete_recent_boosting_scores.bed.gz"),
        afr_ancient=temp(data_path + reference_path + "tmp_AFR_complete_ancient_boosting_scores.bed.gz"),
        eas=temp(data_path + reference_path + "tmp_EAS_complete_boosting_scores.bed.gz"),
        eas_recent=temp(data_path + reference_path + "tmp_EAS_complete_recent_boosting_scores.bed.gz"),
        eas_ancient=temp(data_path + reference_path + "tmp_EAS_complete_ancient_boosting_scores.bed.gz"),
    params:
        url_eur = url_boosting_scores_eur_complete,
        url_eur_recent = url_boosting_scores_eur_recent_complete,
        url_eur_ancient = url_boosting_scores_eur_ancient_complete,
        url_afr= url_boosting_scores_afr_complete,
        url_afr_recent= url_boosting_scores_afr_recent_complete,
        url_afr_ancient= url_boosting_scores_afr_ancient_complete,
        url_eas= url_boosting_scores_eas_complete,
        url_eas_recent=url_boosting_scores_eas_recent_complete,
        url_eas_ancient=url_boosting_scores_eas_ancient_complete
    shell:
        "wget -q -O- {params.url_eur} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.eur}; "
        "wget -q -O- {params.url_eur_recent} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.eur_recent}; "
        "wget -q -O- {params.url_eur_ancient} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.eur_ancient}; "
        "wget -q -O- {params.url_afr} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.afr}; "
        "wget -q -O- {params.url_afr_recent} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.afr_recent}; "
        "wget -q -O- {params.url_afr_ancient} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.afr_ancient}; "
        "wget -q -O- {params.url_eas} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.eas}; "
        "wget -q -O- {params.url_eas_recent} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.eas_recent}; "
        "wget -q -O- {params.url_eas_ancient} | sed 's/\s/\\t/g' | cut -f2-5 | tail -n+2 | gzip > {output.eas_ancient}; "

use rule download_boosting_scores_complete as download_boosting_scores_incomplete with:
    output:
        eur=temp(data_path + reference_path + "tmp_EUR_incomplete_boosting_scores.bed.gz"),
        eur_recent=temp(data_path+ reference_path + "tmp_EUR_incomplete_recent_boosting_scores.bed.gz"),
        eur_ancient=temp(data_path+ reference_path + "tmp_EUR_incomplete_ancient_boosting_scores.bed.gz"),
        afr=temp(data_path + reference_path + "tmp_AFR_incomplete_boosting_scores.bed.gz"),
        afr_recent=temp(data_path + reference_path + "tmp_AFR_incomplete_recent_boosting_scores.bed.gz"),
        afr_ancient=temp(data_path + reference_path + "tmp_AFR_incomplete_ancient_boosting_scores.bed.gz"),
        eas=temp(data_path + reference_path + "tmp_EAS_incomplete_boosting_scores.bed.gz"),
        eas_recent=temp(data_path + reference_path + "tmp_EAS_incomplete_recent_boosting_scores.bed.gz"),
        eas_ancient=temp(data_path + reference_path + "tmp_EAS_incomplete_ancient_boosting_scores.bed.gz"),
    params:
        url_eur = url_boosting_scores_eur_incomplete,
        url_eur_recent = url_boosting_scores_eur_recent_incomplete,
        url_eur_ancient = url_boosting_scores_eur_ancient_incomplete,
        url_afr= url_boosting_scores_afr_incomplete,
        url_afr_recent= url_boosting_scores_afr_recent_incomplete,
        url_afr_ancient= url_boosting_scores_afr_ancient_incomplete,
        url_eas= url_boosting_scores_eas_incomplete,
        url_eas_recent=url_boosting_scores_eas_recent_incomplete,
        url_eas_ancient=url_boosting_scores_eas_ancient_incomplete

use rule lift_1kgp_mask_from_hg19_to_hg38 as lift_boosting_scores_from_hg19_to_hg38 with:
    input:
        bed=data_path + reference_path + "tmp_{population}_{type}_boosting_scores.bed.gz",
        chain=data_path + reference_path + "hg19ToHg38.over.chain.gz"
    output:
        temp(data_path + reference_path + "{population}_{type}_boosting_scores.bed")
    wildcard_constraints:
        pop="EUR|AFR|EAS",
        type="incomplete|complete|incomplete_recent|incomplete_ancient|complete_recent|complete_ancient"

use rule compress_and_index_bstat as compress_and_index_boosting_scores with:
    input:
        data_path + reference_path + "{population}_{type}_boosting_scores.bed"
    output:
        bed=data_path + reference_path + "{population}_{type}_boosting_scores.bed.gz",
        tbi=data_path + reference_path + "{population}_{type}_boosting_scores.bed.gz.tbi"
    wildcard_constraints:
        pop="EUR|AFR|EAS",
        type="incomplete|complete|incomplete_recent|incomplete_ancient|complete_recent|complete_ancient"
