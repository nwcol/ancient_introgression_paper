#!/bin/bash
set -e

chrom=$1
vcf_out=$2
stats_out=$3

u_bar="1.15252e-08"

python simulate_control.py \
    $chrom \
    sexavg_chr${chrom}.txt.gz \
    roulette_10kb_chr${chrom}.csv.gz \
    mask_chr${chrom}.bed.gz \
    intervals_chr${chrom}.txt \
    graph.yaml \
    PopX,PopY \
    $vcf_out \
    $stats_out
