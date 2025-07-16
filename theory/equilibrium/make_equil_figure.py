
import numpy as np
import moments
import msprime
import os.path
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt 


mpl.rcParams["xtick.labelsize"] = 12
mpl.rcParams["ytick.labelsize"] = 12
mpl.rcParams["legend.fontsize"] = 12
mpl.rcParams["font.size"] = 12
mpl.rcParams["axes.titlesize"] = 12
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{amssymb}"
mpl.rcParams["font.family"] = "serif"
mpl.rcParams['font.serif'] = "Computer Modern"
mpl.rcParams["savefig.bbox"] = "tight"


def simulate():
    Ne = 10000
    n = 1000 
    rhos = (1e-2, 1.44, 1e2)
    data = dict()
    for rho in rhos:
        r = rho / (4 * Ne)
        Ts = np.zeros((n, 2))
        for ii in range(n):
            ts = msprime.sim_ancestry(1, discrete_genome=True, 
                recombination_rate=r, sequence_length=2, population_size=Ne)
            tslist = ts.aslist()
            if len(tslist) == 1:
                Ts[ii, :] = tslist[0].get_tmrca(0, 1)
            else:
                Ts[ii, 0] = tslist[0].get_tmrca(0, 1)
                Ts[ii, 1] = tslist[1].get_tmrca(0, 1)
        Ts = Ts / (2 * Ne)
        data[rho] = Ts
    with open("simulation_results.pkl", "wb") as fout:
        pickle.dump(data, fout)
    return 


if not os.path.isfile("simulation_results.pkl"):
    simulate()
with open("simulation_results.pkl", "rb") as fin:
    data = pickle.load(fin)


fig, axes = plt.subplot_mosaic("AAABBB;CCDDEE", figsize=(7, 5.2), 
    layout="constrained", height_ratios=[1.3, 1])


rhos = np.logspace(-8, -1, 50) * 4 * 1e4
y = moments.LD.Demographics1D.snm(rho=rhos, theta=0.001)
Dp = y.H2(0)


axes["A"].plot(rhos, Dp, color="black")
axes["A"].set_ylim(0, 2.1e-6)
axes["A"].set_xscale("log")
axes["A"].set_ylabel(r"$\mathbb{E}[D^+]$")
axes["A"].set_xlabel(r"$\rho$")


CovTT = Dp / 0.001 ** 2 - 1
axes["B"].plot(rhos, CovTT, color="black")
axes["B"].set_ylim(0, 1.05)
axes["B"].set_xscale("log")
axes["B"].set_ylabel("$Cov(T_x, T_y)$")
axes["B"].set_xlabel(r"$\rho$")


axes["C"].scatter(data[1e-2][:, 0], data[1e-2][:, 1], color="black", 
    marker="o", s=4, linewidths=0)
axes["C"].set_ylim(0, 4)
axes["C"].set_xlim(0, 4)
axes["C"].set_xticks([0, 4])
axes["C"].set_yticks([0, 4])
axes["C"].set_ylabel("$T_y$")
axes["C"].set_xlabel("$T_x$")
axes["C"].set_title(r"$\rho=0.01$")


axes["D"].scatter(data[1.44][:, 0], data[1.44][:, 1], color="black", 
    marker="o", s=4, linewidths=0)
axes["D"].set_ylim(0, 4)
axes["D"].set_xlim(0, 4)
axes["D"].set_xticks([0, 4])
axes["D"].set_yticks([0, 4])
axes["D"].set_xlabel("$T_x$")
axes["D"].set_title(r"$\rho=1.44$")


axes["E"].scatter(data[1e2][:, 0], data[1e2][:, 1], color="black", 
    marker="o", s=4, linewidths=0)
axes["E"].set_ylim(0, 4)
axes["E"].set_xlim(0, 4)
axes["E"].set_xticks([0, 4])
axes["E"].set_yticks([0, 4])
axes["E"].set_xlabel("$T_x$")
axes["E"].set_title(r"$\rho=100$")


for label, ax in axes.items():
    ax.set_title(label, loc='left', fontsize='large')


plt.savefig("figure_theory_equilibrium.png", dpi=244)
