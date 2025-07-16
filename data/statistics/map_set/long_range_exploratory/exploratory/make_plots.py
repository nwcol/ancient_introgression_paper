
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


fig, ax = plt.subplots(figsize=(6, 5), layout="constrained")
ax.plot(pairs.sum(0) / pairs.sum(), marker="x")
ax.set_ylabel(r"% pairs")
ax.set_xlabel("bin")
ax.set_yscale("log")
ax.set_xlim(0, pairs.shape[1])
plt.savefig("fig_0_bin_contents.png", dpi=244)


fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.plot(pairs[:, -1] / pairs.sum(1), marker="x")
ax.set_ylabel(r"% of window locus pairs in bin -1")
ax.set_xlabel("chromosome, window")
ax.set_xticks(chrom_start_idxs, chroms)
ax.set_ylim(0,)
ax.set_xlim(0, chrom_start_idxs[-1] + 1)
plt.savefig("fig_1_window_fracs.png", dpi=244)


fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.plot(pairs[:, -1] / pairs[:, -1].sum(), marker="x")
ax.set_ylabel(r"% of bin -1 locus pairs in window")
ax.set_xlabel("chromosome, window")
ax.set_xticks(chrom_start_idxs, chroms)
ax.set_ylim(0,)
ax.set_xlim(0, chrom_start_idxs[-1] + 1)
plt.savefig("fig_2.1_bin-1_fracs.png", dpi=244)



fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.plot(pairs[:, -1] / pairs.sum(), marker="x")
ax.set_ylabel(r"% of all locus pairs in bin -1")
ax.set_xlabel("chromosome, window")
ax.set_xticks(chrom_start_idxs, chroms)
ax.set_ylim(0,)
ax.set_xlim(0, chrom_start_idxs[-1] + 1)
plt.savefig("fig_2_total_fracs.png", dpi=244)


Hs = (sums[:, -1] / sites[:, None])[:, -1]

fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.plot(Hs, marker="x")
ax.set_ylabel(r"window $H$")
ax.set_xlabel(r"chromosome, window")
ax.set_xticks(chrom_start_idxs, chroms)
ax.set_ylim(0,)
ax.set_xlim(0, chrom_start_idxs[-1] + 1)
plt.savefig("fig_3_window_het.png", dpi=244)

fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.scatter(pairs[:, -1] / pairs[:, -1].sum(), Hs ** 2, marker="x")
ax.set_ylabel(r"window $H^2$")
ax.set_xlabel(r"% of window locus pairs in bin -1")
ax.set_ylim(0,)
plt.savefig("fig_4_het_vs_content.png", dpi=244)


Dp_tot = sums[:, :-1, -1].sum(1) / pairs.sum(1)

fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.scatter(pairs[:, -1] / pairs[:, -1].sum(), Dp_tot, marker="x")
ax.set_ylabel(r"window $D^+_{tot}$")
ax.set_xlabel(r"% of window locus pairs in bin -1")
ax.set_ylim(0,)
plt.savefig("fig_5_D+tot_vs_content.png", dpi=244)

uu = mut_facs[:, :-1].sum(1) / pairs.sum(1)

fig, ax = plt.subplots(figsize=(10, 5), layout="constrained")
ax.plot(uu, marker="x")
ax.set_ylabel(r"window $u_Lu_R$")
ax.set_xlabel("chromosome, window")
ax.set_xticks(chrom_start_idxs, chroms)
ax.set_xlim(0, chrom_start_idxs[-1] + 1)
plt.savefig("fig_6_ulur_vs_content.png", dpi=244)