#!/bin/bash
set -e

# naive 
for rep in {0..99}
do
    python naive_bootstrap.py \
        stats/naive_rep_${rep}_stats.pkl \
        raw_stats/empirical_method_0_rep_${rep}
done

# method 0 (adjusts "denominator")
for rep in {0..99}
do
    python weighted_bootstrap_method_0.py \
        stats/method_0_rep_${rep}_stats.pkl \
        raw_stats/empirical_method_0_rep_${rep}
done

# method 1 (adjusts "numerator")
for rep in {0..99}
do
    python naive_bootstrap.py \
        stats/method_1_rep_${rep}_stats.pkl \
        raw_stats/empirical_method_1_rep_${rep}
done

