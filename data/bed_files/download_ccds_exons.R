library(biomaRt)
library(GenomicRanges)
library(tidyverse)
library(data.table)

archives <- listEnsemblArchives()
ensembl_human <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                         dataset="hsapiens_gene_ensembl",
                         host="https://grch37.ensembl.org")
ccds_genes <- getBM(mart=ensembl_human,
                    filters=c("chromosome_name", "with_ccds"),
                    values=list(c(1:22), T),
                    attributes=c("ensembl_gene_id",
                                 "ensembl_exon_id",
                                 "ensembl_transcript_id",
                                 "chromosome_name",
                                 "transcript_start",
                                 "transcript_end",
                                 "exon_chrom_start",
                                 "exon_chrom_end",
                                 "5_utr_start",
                                 "5_utr_end",
                                 "3_utr_start",
                                 "3_utr_end",
                                 "cds_start",
                                 "cds_end",
                                 "cds_length",
                                 "strand"))

fwrite(ccds_genes, "ccds_exons/ccds_genes_hg37.csv")
ccds_genes <- fread("ccds_exons/ccds_genes_hg37.csv")
used_genes <- ccds_genes # ccds (VEP assigns too many MODIFIER sites otherwise)

# the goal is to have a BED file with positions that are CDS in SOME transcript
exon_singletons <- used_genes[!duplicated(used_genes$ensembl_exon_id),]

length(unique(exon_singletons$ensembl_gene_id))
length(unique(used_genes$ensembl_gene_id))

# focus on exons on the tail(s) of each gene
tails <- exon_singletons[!is.na(exon_singletons$`5_utr_start`) |
                           !is.na(exon_singletons$`3_utr_start`),]
tails_pos_strand <- filter(tails, strand==1) # positive strand 5'->3'
tails_neg_strand <- filter(tails, strand==-1) # negative strand 3'->5'

# clipping exon starts
tails_pos_strand[!is.na(tails_pos_strand$`5_utr_end`),]$exon_chrom_start <-
  tails_pos_strand[!is.na(tails_pos_strand$`5_utr_end`),]$`5_utr_end`
tails_pos_strand[!is.na(tails_pos_strand$`3_utr_start`),]$exon_chrom_end <-
  tails_pos_strand[!is.na(tails_pos_strand$`3_utr_start`),]$`3_utr_start`

tails_neg_strand[!is.na(tails_neg_strand$`3_utr_end`),]$exon_chrom_start <-
  tails_neg_strand[!is.na(tails_neg_strand$`3_utr_end`),]$`3_utr_end`
tails_neg_strand[!is.na(tails_neg_strand$`5_utr_start`),]$exon_chrom_end <-
  tails_neg_strand[!is.na(tails_neg_strand$`5_utr_start`),]$`5_utr_start`

# putting it back together
new_tails <- rbind.data.frame(tails_pos_strand, tails_neg_strand)
length(unique(new_tails$ensembl_exon_id)) == nrow(new_tails)

inner_exons <- filter(used_genes, !(ensembl_exon_id %in% new_tails$ensembl_exon_id))
unique_inner_exons <- filter(inner_exons, !duplicated(ensembl_exon_id))
clipped_genes <- rbind.data.frame(new_tails, unique_inner_exons)
clipped_genes <- filter(clipped_genes, exon_chrom_start!=exon_chrom_end)

nrow(clipped_genes) == length(unique(clipped_genes$ensembl_exon_id))
length(unique(clipped_genes$ensembl_gene_id))

clipped_genes <- dplyr::arrange(.data=clipped_genes,
                                chromosome_name,
                                ensembl_gene_id,
                                ensembl_transcript_id,
                                exon_chrom_start,
                                exon_chrom_end)

length(unique(ccds_genes$ensembl_gene_id)) # ccds genes from ensembl
length(unique(clipped_genes$ensembl_gene_id)) # same, after trimming exon ends

for(c in unique(clipped_genes$chromosome_name)) {
  chr_bed <- dplyr::select(clipped_genes, c(chromosome_name,
                                             exon_chrom_start,
                                             exon_chrom_end)) %>%
      dplyr::filter(., chromosome_name == c)
  gr <- makeGRangesFromDataFrame(chr_bed, keep.extra.columns=T) # GR is 1-based
  chr_bed <- as.data.table(GenomicRanges::reduce(gr)) # remove overlaps
  chr_bed <- dplyr::select(chr_bed, c(seqnames, start, end))
  names(chr_bed) <- c("#chrom", "start", "end")

  # to 0-based, half-open intervals
  chr_bed$start <- chr_bed$start - 1
  fwrite(chr_bed, paste("ccds_exons/exons_chr", c, ".bed", sep=""), sep="\t", col.names=F)
}