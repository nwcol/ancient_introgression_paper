
import pickle
import numpy as np
import matplotlib.pyplot as plt


fnames = [f"../raw_stats/Behrer_stats_chr{c}.pkl" for c in range(1, 23)]
regions = {}
for fname in fnames:
    with open(fname, "rb") as fin:
        chrom_stats = pickle.load(fin)
    for label in chrom_stats:
        regions[label] = chrom_stats[label]
sums = np.array([regions[reg]["sums"] for reg in regions])
denoms = np.array([regions[reg]["denoms"] for reg in regions])
sites = denoms[:, -1]
pairs = denoms[:, :-1]
mut_facs = np.array([regions[reg]["mut_facs"] for reg in regions])


labels = list(regions.keys())
chroms = []
chrom_start_idxs = []
for ii, label in enumerate(labels):
    if label[0] not in chroms:
        chroms.append(label[0])
        chrom_start_idxs.append(ii)


H = sums[:, -1, -1] / sites
u = mut_facs[:, -1] / sites
Dptot = sums[:, :-1, -1].sum(1) / pairs.sum(1)
uu_tot = mut_facs[:, :-1].sum(1) / pairs.sum(1)
fig, ax = plt.subplots(figsize=(5, 5), layout="constrained")
ax.scatter(u, H, marker="x")
ax.set_xlabel("$u$")
ax.set_ylabel("$H$")
plt.savefig("fig_0_scatter.png", dpi=244)


Dptot = sums[:, :-1, -1].sum(1) / pairs.sum(1)
uu_tot = mut_facs[:, :-1].sum(1) / pairs.sum(1)
fig, ax = plt.subplots(figsize=(5, 5), layout="constrained")
ax.scatter(uu_tot, Dptot, marker="x")
ax.set_xlabel("$u_lu_r$")
ax.set_ylabel("$D^+_{tot}$")
plt.savefig("fig_1_scatter.png", dpi=244)



means = []
arms = []
last_chrom = 1
last_end = 121500001        
running_num = 0
running_denom = 0
for ii, label in enumerate(labels):
    chrom, _, window = label 
    end = window[-1]
    if chrom != last_chrom or end != last_end:
        print(last_chrom, last_end)
        means.append(running_num / running_denom)
        running_num = 0
        running_denom = 0
        arms.append(chrom)
    running_num += regions[label]["sums"][:-1, -1]
    running_denom += regions[label]["denoms"][:-1]
    last_chrom = chrom
    last_end = end
means.append(running_num / running_denom)
arms.append(chrom)


bins = regions[labels[0]]["bins"]
x = bins[:-1] + np.diff(bins) / 2
fig, _axs = plt.subplots(5, 8, figsize=(14, 8), layout="constrained",
    sharex=True, sharey=False)
axs = _axs.flat
for i in range(39):
    ax = axs[i]
    title = str(arms[i])
    ax.plot(x, means[i], color="tab:blue")
    ax.set_xscale("log")
    ax.set_title(title)

ax = axs[-1]
ax.set_xscale("log")
ax.set_title("TOTAL")

plt.savefig("figure_arms.png", dpi=244)
