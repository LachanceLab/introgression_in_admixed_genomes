#!/usr/bin/env Rscript

library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
args = commandArgs(trailingOnly=TRUE)
lookup <- read.csv(args[1])


ids <- getBM(attributes = c('ensembl_transcript_id_version',
                            'ensembl_gene_id',
                            'ensembl_transcript_id',
                            'ensembl_peptide_id', 
                            'entrezgene_id'),
            filters = 'ensembl_transcript_id_version',
            values = lookup,
            mart = mart)

write.table(ids, file=args[2], row.names=FALSE, quote=FALSE, sep="\t")
