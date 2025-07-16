
import pickle
import numpy as np
import random


def compute_H_sqr(regions):
    # this itself is an approximation~~~ the real one uses pair counts I think
    num = 0
    denom = 0 
    last_chrom = None
    last_end = None
    running_num = 0
    running_denom = 0
    for label in regions:
        chrom, _, window = label 
        end = window[-1]
        
        if chrom != last_chrom or end != last_end:
            num += running_num ** 2
            denom += running_denom ** 2
            running_num = 0
            running_denom = 0

        running_num += regions[label]["sums"][-1]
        running_denom += regions[label]["denoms"][-1]

        last_chrom = chrom
        last_end = end

    return num / denom


def weighted_avg(regions):
    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])
    mut_facs = np.array([regions[reg]["mut_facs"] for reg in regions])

    num_pairs = denoms[:, :-1].sum(1)[:, None]
    tot_pairs = denoms[:, :-1].sum(1)
    tot_facs = mut_facs[:, :-1].sum(1)
    chrom_facs = mut_facs[:, :-1] / (tot_facs / tot_pairs)[:, None]

    norm_stats = np.zeros((sums.shape[1], sums.shape[2]))
    norm_stats[:-1] = sums[:, :-1].sum(0) / chrom_facs.sum(0)[:, None]
    # Naive H 
    norm_stats[-1] = sums[:, -1].sum(0) / denoms[:, -1].sum()
    return norm_stats


def weighted_avg1(regions):
    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])
    mut_facs = np.array([regions[reg]["mut_facs"] for reg in regions])

    facs = (mut_facs[:, :-1].sum(0) / (mut_facs[:, :-1].sum() / denoms[:, :-1].sum()))[:, None]

    norm_stats = np.zeros((sums.shape[1], sums.shape[2]))
    norm_stats[:-1] = sums[:, :-1].sum(0) / facs
    return norm_stats


def naive_avg(regions):
    # No weighting
    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])

    norm_stats = np.zeros((sums.shape[1], sums.shape[2]))
    norm_stats[:-1] = sums[:, :-1].sum(0) / denoms[:, :-1, None].sum(0)
    return norm_stats


stats = {}
for chrom in range(1, 23):
    with open(f"../raw_stats/Behrer_stats_chr{chrom}.pkl", "rb") as fin:
        chrom_stats = pickle.load(fin)
    for label in chrom_stats:
        stats[label] = chrom_stats[label]


means = weighted_avg(stats)


labels = list(stats.keys())
bootmeans = []
for ii in range(len(stats)):
    samples = random.choices(labels, k=len(labels))
    sample_set = {label: stats[label] for label in samples}
    bootmeans.append(weighted_avg(sample_set))
bootmeans = np.array(bootmeans)
varcovs = []
for i in range(len(means)):
    bin_means = np.array([_means[i] for _means in bootmeans])
    varcov_matrix = np.cov(bin_means.T)
    varcovs.append(varcov_matrix)


example = stats[next(iter(stats))]
data = {
    "means": [m for m in means],
    "varcovs": varcovs,
    "pop_ids": example["pop_ids"],
    "bins": np.loadtxt(example["bins"])
}
with open("reweighted.pkl", "wb") as fout:
    pickle.dump(data, fout)
