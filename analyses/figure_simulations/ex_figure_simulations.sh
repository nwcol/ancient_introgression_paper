#!/bin/bash
set -e
graph=$1
pops=$2
out=$3
python simulate_windows_h2.py \
    -g $graph \
    -p $pops \
    -o $out \
    -u 1.3e-8 \
    -r 1e-8 \
    -L 1e6 \
    -n 500
rm *vcf
