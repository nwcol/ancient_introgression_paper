"""
Obtain 1000 bootstrap replicates, compute covariance matrices and save them.
"""    

import pickle
import numpy as np

from dpluspy import bootstrapping


fnames = [f"../raw_stats/Behrer_stats_chr{c}.pkl" for c in range(1, 23)]
regions = {}
for fname in fnames:
    with open(fname, "rb") as fin:
        chrom_stats = pickle.load(fin)
    for label in chrom_stats:
        regions[label] = chrom_stats[label]
bins = regions[label]["bins"]
pop_ids = regions[label]["pop_ids"]
replicates = bootstrapping.get_bootstrap_reps(
    regions, num_reps=1000, weighted=True)
varcovs = bootstrapping.compute_varcovs(replicates)
means = bootstrapping.weighted_means_across_regions(regions)

# Subset statistics
subset_reps = []
for rep_means in replicates:
    rep = {"bins": bins, "pop_ids": pop_ids, "means": rep_means, 
           "varcovs": varcovs}
    sub_bins, sub_means, sub_varcovs = bootstrapping.subset_stats(
        rep, min_r=1e-6, max_r=1e-2, return_dict=False)
    subset_reps.append(sub_means)

overall = {"bins": bins, "pop_ids": pop_ids, "means": means, 
        "varcovs": varcovs}
sub_bins, sub_means, sub_varcovs = bootstrapping.subset_stats(
    overall, min_r=1e-6, max_r=1e-2, return_dict=False)

data = {
    "means": sub_means,
    "varcovs": sub_varcovs,
    "replicates": subset_reps,
    "bins": sub_bins,
    "pop_ids": pop_ids
}
with open("Bootstrap_BehrerStats16bins.pkl", "wb") as fout:
    pickle.dump(data, fout)

