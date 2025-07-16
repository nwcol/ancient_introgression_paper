"""
Plot chromosome arm/total means and bootstrap CIs of these statistics:
H^2 (naive), H^2 (properly weighted), D+tot (naive), D+long (naive),
D+tot (weighted), D+long (naive).

With our window-overhang architecture, there *is no* way to adequately bootstrap
the `correctly` weighted H^2 statistic!

We just plot the Yoruba3 statistics for now for simplicity.
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
import random


def load_stats(idx):
    fnames = [f"../raw_stats/Behrer_stats_chr{c}.pkl" for c in range(1, 23)]
    regions = {}
    for fname in fnames:
        with open(fname, "rb") as fin:
            chrom_stats = pickle.load(fin)
        for label in chrom_stats:
            regions[label] = chrom_stats[label]
    for label in regions:
        regions[label]["sums"] = regions[label]["sums"][:, [idx]].flatten()
    return regions


def separate_chrom_arms(regions):
    labels = list(regions.keys())
    arm_sets = dict()
    last_chrom = 0
    last_end = 0
    for ii, label in enumerate(labels):
        chrom, _, window = label 
        end = window[-1]
        if chrom != last_chrom:
            key = (chrom, 0)
            arm_sets[key] = dict()
        elif end != last_end:
            key = (chrom, 1)
            arm_sets[key] = dict()
        arm_sets[key][label] = regions[label]
        last_chrom = chrom
        last_end = end
    return arm_sets


def get_arm_avg(regions):
    # compute: H^2naive, H^2, D+longnaive, D+totnaive, D+longw, D+totw.
    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])
    mut_facs = np.array([regions[reg]["mut_facs"] for reg in regions])

    tot_pairs = denoms[:, :-1].sum(1)
    tot_facs = mut_facs[:, :-1].sum(1)
    chrom_facs = mut_facs[:, :-1] / (tot_facs / tot_pairs)[:, None]

    stats = np.zeros(6)

    stats[0] = (sums[:, -1].sum() / denoms[:, -1].sum()) ** 2
    stats[1] = stats[1]
    stats[2] = sums[:, :-1].sum() / denoms[:, :-1].sum()
    if denoms[:, -2].sum() > 0:
        stats[3] = sums[:, -2].sum() / denoms[:, -2].sum()
    else:
        idx = np.where(denoms[:, :-1].sum(0) > 0)[0][-1]
        stats[3] = sums[:, idx].sum() / denoms[:, idx].sum()
    stats[4] = sums[:, :-1].sum() / chrom_facs.sum()
    if denoms[:, -2].sum() > 0:
        stats[5] = sums[:, -2].sum() / chrom_facs[:, -1].sum()
    else:
        stats[5] = sums[:, idx].sum() / chrom_facs[:, idx].sum()

    return stats


def get_tot_avg(regions, arm_stats):

    sums = np.array([regions[reg]["sums"] for reg in regions])
    denoms = np.array([regions[reg]["denoms"] for reg in regions])
    mut_facs = np.array([regions[reg]["mut_facs"] for reg in regions])

    tot_pairs = denoms[:, :-1].sum(1)
    tot_facs = mut_facs[:, :-1].sum(1)
    chrom_facs = mut_facs[:, :-1] / (tot_facs / tot_pairs)[:, None]

    stats = np.zeros(6)

    stats[2] = sums[:, :-1].sum() / denoms[:, :-1].sum()
    stats[3] = sums[:, -2].sum() / denoms[:, -2].sum()
    stats[4] = sums[:, :-1].sum() / chrom_facs.sum()
    stats[5] = sums[:, -2].sum() / chrom_facs[:, -1].sum()

    stats[0] = (sums[:, -1].sum() / denoms[:, -1].sum()) ** 2

    num = 0
    denom = 0
    for key in arm_stats:
        num += np.array([arm_stats[key][x]["sums"] for x in arm_stats[key]])[:, -1].sum() ** 2
        denom += np.array([arm_stats[key][x]["denoms"] for x in arm_stats[key]])[:, -1].sum() ** 2
    stats[1] = num / denom

    return stats


def run_bootstrap(regions, arm_sums=None):
    if arm_sums is None:
        means = get_arm_avg(regions)
    else:
        means = get_tot_avg(regions, arm_sums)
    labels = list(regions.keys())
    bootmeans = []
    for ii in range(len(regions)):
        samples = random.choices(labels, k=len(labels))
        sample_set = {label: regions[label] for label in samples}
        sample_means = get_arm_avg(sample_set)
        bootmeans.append(sample_means)
    bootmeans = np.array(bootmeans)
    stddevs = []
    for i in range(len(means)):
        bin_means = np.array([_means[i] for _means in bootmeans])
        bin_stds = np.std(bin_means, axis=0)
        stddevs.append(bin_stds)
    return means, np.array(stddevs)


regions = load_stats(-1)
# Compute arm statistics
arm_sums = separate_chrom_arms(regions)
labels = list(arm_sums.keys())
arm_stats = {}
for key in arm_sums:
    arm_stats[key] = run_bootstrap(arm_sums[key])
arm_means = np.array([arm_stats[x][0] for x in arm_stats])
arm_errs = 1.96 * np.array([arm_stats[x][1] for x in arm_stats])
tot_means, tot_stds = run_bootstrap(regions, arm_sums=arm_sums)
tot_errs = 1.96 * tot_stds
    # compute: H^2naive, H^2, D+longnaive, D+totnaive, D+longw, D+totw.

fig, _axs = plt.subplots(5, 8, figsize=(14, 8), layout="constrained",
    sharex=True, sharey=False)
axs = _axs.flat

for i, key in enumerate(arm_stats):
    ax = axs[i]
    title = f'{labels[i][0]}_{labels[i][1]}'
    ax.set_title(title)

    ax.errorbar([0], arm_means[i, 0], arm_errs[i, 0], capsize=2, color="tab:blue", marker="x")
    ax.errorbar([1], arm_means[i, 2], arm_errs[i, 2], capsize=2, color="tab:orange", marker="o")
    ax.errorbar([1.2], arm_means[i, 3], arm_errs[i, 3], capsize=2, color="tab:orange", marker="o", markerfacecolor="none",)
    ax.errorbar([2], arm_means[i, 4], arm_errs[i, 4], capsize=2, color="tab:green", marker="o")
    ax.errorbar([2.2], arm_means[i, 5], arm_errs[i, 5], capsize=2, color="tab:green", marker="o", markerfacecolor="none",)
    ax.grid(alpha=0.3)

ax = axs[-1]
ax.set_title("Total")
l1 = ax.errorbar([0], tot_means[0], tot_errs[0], capsize=2, color="tab:blue", marker="x", label="naive $H^2$")
l2 = ax.errorbar([0.2], tot_means[1], tot_errs[1], capsize=2, color="tab:red", marker="x", label="$H^2$")
l3 = ax.errorbar([1], tot_means[2], tot_errs[2], capsize=2, color="tab:orange", marker="o", label="naive $D^+_{tot}$")
l4 = ax.errorbar([1.2], tot_means[3], tot_errs[3], capsize=2, color="tab:orange", marker="o", markerfacecolor="none", label="naive $D^+_{long}$")
l5 = ax.errorbar([2], tot_means[4], tot_errs[4], capsize=2, color="tab:green", marker="o", label="$D^+_{tot}$")
l6 = ax.errorbar([2.2], tot_means[5], tot_errs[5], capsize=2, color="tab:green", marker="o", markerfacecolor="none", label="$D^+_{long}$")
ax.set_xticks([], [])
ax.grid(alpha=0.3)

fig.legend(framealpha=0, ncols=6, loc="lower center", bbox_to_anchor=(0.5, -0.1))
plt.savefig("fig_arm_statistics", dpi=244, bbox_inches="tight")

















