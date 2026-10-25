#!/bin/bash
set -e
unzip data.zip
python simulate_genome.py $@
rm data_chr*.vcf

