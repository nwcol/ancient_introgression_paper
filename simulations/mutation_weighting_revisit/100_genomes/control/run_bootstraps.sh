#!/bin/bash
set -e

for rep in {0..99}
do
    python bootstrap.py \
        stats/control_rep_${rep}_stats.pkl \
        raw_stats/control_rep_${rep}
done
