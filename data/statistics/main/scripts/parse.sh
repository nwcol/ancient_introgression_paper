#!/bin/bash
set -e
chrom=$1
gzip -d "roulette_b37_chr${chrom}.npy.gz"
python parse.py $chrom
rm "roulette_b37_chr${chrom}.npy"