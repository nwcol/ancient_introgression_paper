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
pred_genes <- getBM(mart=ensembl_human,
                    filters=c("chromosome_name"),
                    values=list(1:22),
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
fwrite(ccds_genes, "ccds_genes_hg38.csv")
fwrite(pred_genes, "pred_genes_hg38.csv")
ccds_genes <- fread("ccds_genes_hg38.csv")
pred_genes <- fread("pred_genes_hg38.csv")
# the goal is to have a BED file with positions that are CDS in SOME transcript
exon_singletons <- pred_genes[!duplicated(pred_genes$ensembl_exon_id),]
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
inner_exons <- filter(pred_genes,
                      !(ensembl_exon_id %in% new_tails$ensembl_exon_id))
unique_inner_exons <- filter(inner_exons, !duplicated(ensembl_exon_id))
clipped_genes <- rbind.data.frame(new_tails, unique_inner_exons)
clipped_genes <- clipped_genes[clipped_genes$exon_chrom_end !=
                               clipped_genes$exon_chrom_start,]
nrow(clipped_genes) == length(unique(clipped_genes$ensembl_exon_id))
clipped_genes <- dplyr::arrange(.data=clipped_genes,
                                chromosome_name,
                                ensembl_gene_id,
                                ensembl_transcript_id,
                                exon_chrom_start,
                                exon_chrom_end)
clipped_genes <- filter(clipped_genes, exon_chrom_start!=exon_chrom_end)
# this discards most genes that are not CCDS
genes_sh <- merge(clipped_genes, sh_table, by="ensembl_gene_id")
genes_sh <- dplyr::arrange(.data=genes_sh, chromosome_name, exon_chrom_start)
sum(genes_sh$exon_chrom_end - genes_sh$exon_chrom_start)
# check if we are losing genes when mergin with the Zeng et al table
length(unique(pred_genes$ensembl_gene_id)) # ccds genes from ensembl
length(unique(clipped_genes$ensembl_gene_id)) # same, after trimming exon ends
length(unique(sh_table$ensembl_gene_id)) # zeng et al table for autosomes
length(unique(genes_sh$ensembl_gene_id)) # merged table
clipped_leftover <- filter(clipped_genes,
                           !ensembl_gene_id %in% genes_sh$ensembl_gene_id)
ccds_clipped_leftover <- filter(clipped_leftover,
                                ensembl_gene_id %in% ccds_genes$ensembl_gene_id)
length(unique(ccds_clipped_leftover$ensembl_gene_id))
# including genes not in Zeng table as 11th decile, but only those with CCDS
ccds_clipped_leftover$hgnc <- NA
ccds_clipped_leftover$chrom <- NA
ccds_clipped_leftover$obs_lof <- NA
ccds_clipped_leftover$exp_lof <- NA
ccds_clipped_leftover$prior_mean <- NA
ccds_clipped_leftover$post_mean <- NA
ccds_clipped_leftover$post_lower_95 <- NA
ccds_clipped_leftover$post_upper_95 <- NA
ccds_clipped_leftover$decile_constraint <- 11 # 11th decile
genes_sh_ext <- rbind.data.frame(genes_sh, ccds_clipped_leftover)
genes_sh_ext <- dplyr::arrange(.data=genes_sh_ext, chromosome_name, exon_chrom_start)
length(unique(genes_sh_ext$ensembl_gene_id)) # 18770
# genes in the Zeng et al table not in the predicted genes from Ensenbl
zeng_leftover <- filter(sh_table, !ensembl_gene_id %in% genes_sh_ext$ensembl_gene_id)
nrow(zeng_leftover) # 275
fwrite(genes_sh_ext, "genes_sh_deciles_hg38.csv", sep=",")
genes_sh_ext <- fread("genes_sh_deciles_hg38.csv")
cds_sums <- numeric(length(unique(genes_sh_ext$chromosome_name)))
for(c in unique(genes_sh_ext$chromosome_name)) {
  chr_bed <- dplyr::select(genes_sh_ext, c(chromosome_name,
                                           exon_chrom_start,
                                           exon_chrom_end)) %>%
             dplyr::filter(., chromosome_name == c)
  gr <- makeGRangesFromDataFrame(chr_bed, keep.extra.columns=T) # GR is 1-based
  chr_bed <- as.data.table(reduce(gr)) # remove overlaps
  chr_bed <- dplyr::select(chr_bed, c(seqnames, start, end))
  names(chr_bed) <- c("#chrom", "start", "end")
  cds_sums[c] <- sum(chr_bed$end - chr_bed$start)
  # to 0-based, half-open intervals
  chr_bed$start <- chr_bed$start - 1
  fwrite(chr_bed, paste("~/Data/bgs_lmr/human_data/data/annotation_maps/bedfiles/exons/chr",
                        c, "_exons.bed", sep=""), sep="\t", col.names=F)
  for(d in unique(genes_sh_ext$decile_constraint)) {
    tbl <- dplyr::filter(genes_sh_ext, chromosome_name==c, decile_constraint==d)
    df <- dplyr::select(tbl, c("chromosome_name",
                               "exon_chrom_start",
                               "exon_chrom_end"))
    gr <- makeGRangesFromDataFrame(df, keep.extra.columns=T) # GR is 1-based
    df <- as.data.table(reduce(gr)) # remove overlaps
    df <- dplyr::select(df, c(seqnames, start, end))
    names(df) <- c("#chrom", "start", "end")
    # to 0-based, half-open intervals
    df$start <- df$start - 1
    fwrite(df, paste("~/Data/bgs_lmr/human_data/data/sfs_genes/masks/bed/", c,
                     "_exons_lof_d", d, ".tsv.gz", sep=""), sep="\t")
  }
}