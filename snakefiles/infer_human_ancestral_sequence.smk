import os
from Bio import SeqIO

envvars:
    "PERL5LIB"

# get the human ancestral sequence
rule get_human_ancestral_sequence:
    output:
        temp(expand(data_path + human_ancestral_sequence_path + 'tmp_homo_sapiens_ancestor_{chr}.fa', chr=chromosomes))
    params:
        ensembl = ensembl_path,
        data = data_path,
        ancestral = human_ancestral_sequence_path,
        x = os.environ['PERL5LIB']
    conda:
        "../envs/ensembl.yaml"
    shell:
        "export PERL5LIB={params.x}; "
        "{params.ensembl}ensembl-compara/scripts/ancestral_sequences/get_ancestral_sequence.pl --dir "
        "{params.data}{params.ancestral}"

rule format_human_ancestral_sequence:
    input:
        data_path + human_ancestral_sequence_path + 'tmp_homo_sapiens_ancestor_{chr}.fa'
    output:
        data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_{chr}.fa'
    params:
        chrom="{chr}"
    shell:
        "cp {input} {output}; sed -i 's/>/>chr{params.chrom} /' {output}"


# convert human ancestral sequence to BED file
rule human_ancestral_sequence_to_bed:
    input:
        data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_{chr}.fa'
    output:
        temp(data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_hg38_{chr}.bed')
    run:
        chrom = wildcards.chr
        sequence = list(SeqIO.parse(input[0], 'fasta'))[0].seq
        output_file = open(output[0], 'w')
        for i, bp in enumerate(sequence):
            if bp == '-' or bp == '.' or bp == 'N':
                continue
            output_file.write(f'chr{chrom}\t{i}\t{i + 1}\t{bp}\n')
        output_file.close()

rule compress_and_index_ancestral_states:
    input:
        data_path+ human_ancestral_sequence_path + 'homo_sapiens_ancestor_hg38_{chr}.bed'
    output:
        multiext(data_path + human_ancestral_sequence_path + 'homo_sapiens_ancestor_hg38_{chr}.bed', '.gz', '.gz.tbi')
    conda:
        "../envs/tabix.yaml"
    shell:
        "bgzip {input}; tabix {input}.gz"