#!/bin/bash
set -e
chrom=$1


python run_simulation.py $chrom


python compute_stats.py $chrom "arms_chr${chrom}.txt" "0" "arm_stats_chr${chrom}.pkl"

python compute_stats.py $chrom "intervals_chr${chrom}.txt" "1" "interval_overhang_stats_chr${chrom}.pkl"

python compute_stats.py $chrom "intervals_chr${chrom}.txt" "0" "interval_stats_chr${chrom}.pkl"