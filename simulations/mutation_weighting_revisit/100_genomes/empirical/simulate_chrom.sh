#!/bin/bash
set -e

chrom=$1
vcf_out=$2
output_0=$3
output_1=$4

u_bar="1.1525e-08"

python simulate_chrom.py \
    $chrom \
    sexavg_chr${chrom}.txt.gz \
    roulette_10kb_chr${chrom}.csv.gz \
    mask_chr${chrom}.bed.gz \
    intervals_chr${chrom}.txt \
    graph.yaml \
    PopX,PopY \
    $vcf_out \
    $output_0

python compute_adj_stats.py \
    $vcf_out \
    mask_chr${chrom}.bed.gz \
    sexavg_chr${chrom}.txt.gz \
    roulette_10kb_chr${chrom}.csv.gz \
    intervals_chr${chrom}.txt \
    $chrom \
    $u_bar \
    $output_1

