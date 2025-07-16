#!/bin/bash
set -e
chrom=$1


python compute_stats.py $chrom "intervals_chr${chrom}.txt" "1" "interval_overhang_stats_chr${chrom}.pkl"