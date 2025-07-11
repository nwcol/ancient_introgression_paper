#!/bin/bash
chrom=$1
gzip -d "roulette_b37_chr${chrom}.npy.gz"
python parse_intergenic.py $@
rm "roulette_b37_chr${chrom}.npy"