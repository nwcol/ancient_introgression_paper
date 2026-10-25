#!/bin/bash
set -e
unzip data_chr$1
python simulate_chromosome.py $@
rm temp.vcf
rm -r data_chr$1
