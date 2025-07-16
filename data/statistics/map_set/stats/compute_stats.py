"""
Perform bootstraps on the 11 data sets. 
"""


import pickle
import numpy as np
import random

from dpluspy import bootstrapping


def weighted_avg(regions):
    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])
    mut_facs = np.array([regions[reg]["mut_facs"] for reg in regions])

    tot_pairs = denoms[:, :-1].sum(1)
    tot_facs = mut_facs[:, :-1].sum(1)
    chrom_facs = mut_facs[:, :-1] / (tot_facs / tot_pairs)[:, None]

    norm_stats = np.zeros((sums.shape[1], sums.shape[2]))
    norm_stats[:-1] = sums[:, :-1].sum(0) / chrom_facs.sum(0)[:, None]
    # Naive H 
    norm_stats[-1] = sums[:, -1].sum(0) / denoms[:, -1].sum()
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


######





def call_bootstrap(fnames, out_fname, weighted=False):
    regions = {}
    for fname in fnames:
        with open(fname, "rb") as fin:
            chrom_stats = pickle.load(fin)
        for label in chrom_stats:
            regions[label] = chrom_stats[label]
    example = regions[label]
    means, varcovs = bootstrapping.bootstrap_stats(regions, weighted=weighted)
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
        'hapmap',
        'Hinch',
        'omniFIN',
        'omniLWK',
        'omniYRI',
        'pyrhoFIN',
        'pyrhoLWK',
        'pyrhoYRI',
        'ZhouFHS',
        'ZhouJHS'
    ]
    for map_type in map_types:
        if map_type == "hapmap":
            continue
        fnames = [f"../raw_stats/{map_type}_stats_chr{c}.pkl" for c in range(1, 23)]
        out_fname = f"{map_type}_stats.pkl"
        call_bootstrap(fnames, out_fname, weighted=True)
        out_fname = f"{map_type}_naive_stats.pkl"
        call_bootstrap(fnames, out_fname, weighted=False)


    # hapmap is missing chromosomes 6, 7
    subchroms = list(range(1, 23))
    subchroms.remove(6)
    subchroms.remove(7)
    for map_type in map_types:
        fnames = [f"../raw_stats/{map_type}_stats_chr{c}.pkl" for c in subchroms]
        out_fname = f"{map_type}_skip6_7_stats.pkl"
        call_bootstrap(fnames, out_fname)

    return


main()

