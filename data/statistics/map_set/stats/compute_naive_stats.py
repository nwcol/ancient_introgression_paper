"""
Perform bootstraps on 4 data sets w/o weighting.
"""


import pickle
import numpy as np
import random


def weighted_avg(regions):
    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])
    norm_stats = sums.sum(0) / denoms.sum(0)[:, None]
    return norm_stats


def run_bootstrap(fnames, out_fname):
    regions = {}
    for fname in fnames:
        with open(fname, "rb") as fin:
            chrom_stats = pickle.load(fin)
        for label in chrom_stats:
            regions[label] = chrom_stats[label]
    means = weighted_avg(regions)
    labels = list(regions.keys())
    bootmeans = []
    for ii in range(len(regions)):
        samples = random.choices(labels, k=len(labels))
        sample_set = {label: regions[label] for label in samples}
        bootmeans.append(weighted_avg(sample_set))
    bootmeans = np.array(bootmeans)
    varcovs = []
    for i in range(len(means)):
        bin_means = np.array([_means[i] for _means in bootmeans])
        varcov_matrix = np.cov(bin_means.T)
        varcovs.append(varcov_matrix)
    example = regions[next(iter(regions))]
    data = {
        "means": [m for m in means],
        "varcovs": varcovs,
        "pop_ids": example["pop_ids"],
        "bins": example["bins"]
    }
    with open(out_fname, "wb") as fout:
        pickle.dump(data, fout)
    return


def main():
    map_types = [
        'Behrer',
        'omniYRI',
        'ZhouFHS',
        'ZhouJHS'
    ]
    for map_type in map_types:
        if map_type == "hapmap":
            continue
        fnames = [f"../raw_stats/{map_type}_stats_chr{c}.pkl" for c in range(1, 23)]
        out_fname = f"{map_type}_naive_stats.pkl"
        run_bootstrap(fnames, out_fname)
    return


main()

