#!/bin/bash
set -e
# Grab arguments
g_truth=$1
g_0=$2
p_0=$3
g_1=$4
p_1=$5
prefix=$6
# Construct output names
sim_data="${prefix}_data.pkl"
# Simulate dataset
python simulate_windows_sfs.py \
    --graph $g_truth \
    --pops pops.txt \
    --out $sim_data \
    -u 1.3e-8 \
    -r 1e-8 \
    -L 1e6 \
    -n 500 \
    --n-samples 1
# Delete VCF files
rm *vcf
# Fit nested models serially
for i in {0..2}
do
    sub_prefix="${prefix}_try_${i}"
    python fit_nested_models_sfs.py \
        --data $sim_data \
        --graph0 $g_0 \
        --graph1 $g_1 \
        --params0 $p_0 \
        --params1 $p_1 \
        --nested Tpulse 10000 p12 1e-3 \
        -o $sub_prefix \
        -u 1.3e-8 \
        --perturb 0.30
done
