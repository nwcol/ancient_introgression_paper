"""
Simulate a chromosome with empirical recombination/mutation maps under a given 
demographic model and compute D+ statistics across it.
"""

import demes
import msprime 
import numpy as np
import dpluspy
import pandas
import sys
import pickle


# Arguments; filepaths
chrom = sys.argv[1]
rec_map_fname = sys.argv[2]
mut_map_fname = sys.argv[3]
mask_fname = sys.argv[4]
window_fname = sys.argv[5]
graph_fname = sys.argv[6]
pop_ids = sys.argv[7].split(",")


output_vcf = sys.argv[8]
output_stats = sys.argv[9]


# Load mutation map (a .csv file with uniform windows)
mut_df = pandas.read_csv(mut_map_fname)
edges = np.append(
    np.array(mut_df["chromStart"]), np.array(mut_df["chromEnd"])[-1])
vals = np.array(mut_df["mut_map"])
vals[np.isnan(vals)] = 0
mut_map = msprime.RateMap(position=edges, rate=vals)


# Load recombination map
sequence_length = edges[-1]
rec_map = msprime.RateMap.read_hapmap(
    rec_map_fname, sequence_length=sequence_length, position_col=1, rate_col=2)


# Define bins in `r`
bins = np.logspace(-6, -2, 17)


# Define demography
samples = {pop_id: 1 for pop_id in pop_ids}
graph = demes.load(graph_fname)
demog = msprime.Demography.from_demes(graph)


# Run the coalescent simulation
ts = msprime.sim_ancestry(
    demography=demog,
    samples=samples,
    sequence_length=sequence_length,
    recombination_rate=rec_map,
    record_provenance=False,
    discrete_genome=True
)


# Simulate mutations
mts = msprime.sim_mutations(ts, rate=mut_map, record_provenance=False)
with open(output_vcf, "w") as fout:
    mts.write_vcf(fout, position_transform=dpluspy.utils._increment1)


# Compute statistics
stats = dpluspy.parsing.parse_stats(
    output_vcf,
    bed_file=mask_fname,
    rec_map_file=rec_map_fname,
    pos_col="pos",
    map_col="cM",
    mut_map_file=mut_map_fname,
    mut_col="mut_map",
    interval_file=window_fname,
    r_bins=bins,
    chrom=chrom,
    get_denoms=False
)

with open(output_stats, "wb") as fout:
    pickle.dump(stats, fout)




