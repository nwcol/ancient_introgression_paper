#!/bin/bash
chrom=$1
gzip -d "roulette_b37_chr${chrom}.npy.gz"
python parse_with_map.py $@
rm "roulette_b37_chr${chrom}.npy"