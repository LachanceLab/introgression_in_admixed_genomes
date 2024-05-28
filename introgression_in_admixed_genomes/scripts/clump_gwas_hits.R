#!/usr/bin/env Rscript
library('LDlinkR')
library('optparse')

option_list <- list(
  make_option(c("-c", "--fn_catalog"), type="character", default=NULL,
              help="File path to the catalog file"),
  make_option(c("-s", "--fn_segments"), type="character", default=NULL,
              help="File path to the segments file"),
  make_option(c("-a", "--association_col"), type="character", default=NULL,
              help="Name of the association column"),
  make_option(c("-t", "--token"), type="character", default=NULL,
              help="Token for accessing LDlink"),
  make_option(c("-o", "--fn_out"), type="character", default=NULL,
              help="Output file path"),
  make_option(c("-p", "--fn_overlap"), type="character", default=NULL,
              help="File path to the overlap file")
)
parser <- OptionParser(option_list=option_list)
opt <- parse_args(parser)

header_reference <- names(as.data.frame(read.csv(opt$fn_catalog, header=TRUE, sep='\t')))
header_segments <- names(as.data.frame(read.csv(opt$fn_segments, header=TRUE, sep='\t')))
header_segments[1] <- 'chrom_r'
header_segments[2] <- 'start_r'
header_segments[3] <- 'end_r'
columns <-  c(header_reference, header_segments, 'overlap')
df <- as.data.frame(read.csv(opt$fn_overlap, header=F, col.names=columns, sep='\t'))
clumped <- list()
for (trait in unique(df[, opt$association_col])) {
    c_df <- df[df[opt$association_col] == trait, ]
    c_df <- c_df[!duplicated(c_df["SNPS"]), ]
    if (dim(c_df)[1] > 1) {
        for (chrom in unique(c_df$CHR_ID)) {
            c_chrom <- c_df[c_df$CHR_ID == chrom, ]
            if (dim(c_chrom)[1] > 1) {
                snps <- SNPclip(c_chrom[, "SNPS"], token=opt$token, genome_build="grch38")
                snps <- as.data.frame(snps)
                o_snps <- snps
                snps <- snps[snps$Details == 'Variant kept.', ]
                if (is.null(dim(snps))) {
                    print(o_snps)
                    c_chrom_clumped <- c_chrom[!duplicated(c_chrom[c("chrom_r", 'start_r', 'end_r')]), ]
                }else if (dim(snps)[1] == 1) {
                    c_chrom_clumped <- c_chrom[c_chrom["SNPS"] == snps[, "RS_Number"],]
                }else {
                    c_chrom_clumped <- c_chrom[c_chrom[, "SNPS"] %in% snps$RS_Number, ]
                }
                clumped <- rbind(clumped, c_chrom_clumped)
            }else {
                clumped <- rbind(clumped, c_chrom)
            }
        }
    }else {
        clumped <- rbind(clumped, c_df)
    }
}
clumped <- as.data.frame(clumped)
write.table(clumped, opt$fn_out, row.names=F, col.names=T, quote=F, sep="\t")
