
import demes
import msprime
import numpy as np
import pandas
import pickle
import sys

import h2py
from h2py.utils import timestamp


# Input arguments
chrom = sys.argv[1]
length = int(sys.argv[2])
graph_file = sys.argv[3]
pop_file = sys.argv[4]
out_stem = sys.argv[5]


# Filenames
mut_map_file = f"data_chr{chrom}/roulette_1kb_chr{chrom}.tsv.gz"
bed_file = f"data_chr{chrom}/roulette_coverage_chr{chrom}.bed.gz"
rec_map_file = f"data_chr{chrom}/sexavg_chr{chrom}.txt.gz"


# Setup
graph = demes.load(graph_file)
ws = 5_000_000
intervals = [[i*ws, (i+1)*ws] for i in range(int(length//ws) + 1)]
r_bins = np.logspace(-6, -2, 17)


def run_simulation():
    pop_spec = h2py.matrices._read_pop_file(pop_file)
    samples = {pop: 1 for pop in pop_spec}
    demog = msprime.Demography.from_demes(graph)

    rec_map = msprime.RateMap.read_hapmap(rec_map_file, sequence_length=length)
    mut_tbl = pandas.read_table(mut_map_file, engine="python", sep=None)
    mut_rate = mut_tbl["mutRate"]
    mut_rate[np.isnan(mut_rate)] = 0
    starts = np.array(mut_tbl["chromStart"])
    mut_pos = np.append(starts[starts < length], length)
    mut_map = msprime.RateMap(position=mut_pos, rate=mut_rate[starts < length])

    ts = msprime.sim_ancestry(
        samples,
        demography=demog,
        sequence_length=length,
        recombination_rate=rec_map,
    )
    ts = msprime.sim_mutations(ts, rate=mut_map)
    with open("temp.vcf", "w+") as f:
        ts.write_vcf(f, position_transform = lambda x: 1 + np.array(x))
    return


def compute_stats():
    sums = {}
    for i, interval in enumerate(intervals):
        sums[f"{out_stem}_{i}"] = h2py.parsing.compute_h2_stats(
            vcf_file="temp.vcf",
            pop_file=pop_file,
            interval=interval,
            bed_file=bed_file,
            rec_map_file=rec_map_file,
            rec_pos_col="pos",
            rec_map_col="cM",
            mut_map_file=mut_map_file,
            r_bins=r_bins,
            report=True,
        )
    return sums


def compute_stats_naive():
    sums = {}
    for i, interval in enumerate(intervals):
        sums[f"{out_stem}_{i}"] = h2py.parsing.compute_h2_stats(
            vcf_file="temp.vcf",
            pop_file=pop_file,
            interval=interval,
            bed_file=bed_file,
            rec_map_file=rec_map_file,
            r_bins=r_bins,
            report=True
        )
    return sums


run_simulation()

adj_stats = compute_stats()
adj_output = f"{out_stem}_adj.pkl"
with open(adj_output, "wb") as fout:
    pickle.dump(adj_stats, fout)

naive_stats = compute_stats_naive()
naive_output = f"{out_stem}_naive.pkl"
with open(naive_output, "wb") as fout:
    pickle.dump(naive_stats, fout)

