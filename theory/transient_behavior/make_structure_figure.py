
import numpy as np
import moments
import msprime
import demes
import demesdraw
import pickle
import os
import matplotlib.pyplot as plt 


# Parameters
u = 1.4e-8
rs = np.logspace(-6, -2, 20)
ts = np.linspace(0, 12000, 30)

# for scatter plots
rrs = [1e-5, 1e-4, 4e-4]


b = demes.Builder()
b.add_deme("A", epochs=[{"end_time": 11000, "start_size": 10000}])
b.add_deme("S1", ancestors=["A"], epochs=[{"end_time": 1000, "start_size": 8000}])
b.add_deme("S2", ancestors=["A"], epochs=[{"end_time": 1000, "start_size": 2000}])
b.add_deme("S",  ancestors=["S1", "S2"], proportions=[0.5, 0.5], start_time=1e3,
           epochs=[{"end_time": 0, "start_size": 10000}])
g_struc = b.resolve()
demog_struc = msprime.Demography.from_demes(g_struc)
debug = demog_struc.debug()
str_rates, _ = debug.coalescence_rate_trajectory(ts, {"S": 2})


# Construct rate-equivalent panmictic model
b = demes.Builder()
epochs = []
for t, rate in zip(ts[::-1], str_rates[::-1]):
    epochs.append({"end_time": t, "start_size": 1 / (2 * rate)})
b.add_deme("S", epochs=epochs)
g_pan = b.resolve()
demog_pan = msprime.Demography.from_demes(g_pan)
debug = demog_pan.debug()
pan_rates, _ = debug.coalescence_rate_trajectory(ts, {"S": 2})


def simulate():
    n = 1000
    data = {"struc": {}, "pan": {}}
    for r in rrs:
        Ts = np.zeros((n, 2))
        for ii in range(n):
            ts = msprime.sim_ancestry({"S": 1}, discrete_genome=True, 
                recombination_rate=r, sequence_length=2, demography=demog_struc)
            tslist = ts.aslist()
            if len(tslist) == 1:
                Ts[ii, :] = tslist[0].get_tmrca(0, 1)
            else:
                Ts[ii, 0] = tslist[0].get_tmrca(0, 1)
                Ts[ii, 1] = tslist[1].get_tmrca(0, 1)
        data["struc"][r] = Ts

        Ts = np.zeros((n, 2))
        for ii in range(n):
            ts = msprime.sim_ancestry({"S": 1}, discrete_genome=True, 
                recombination_rate=r, sequence_length=2, demography=demog_pan)
            tslist = ts.aslist()
            if len(tslist) == 1:
                Ts[ii, :] = tslist[0].get_tmrca(0, 1)
            else:
                Ts[ii, 0] = tslist[0].get_tmrca(0, 1)
                Ts[ii, 1] = tslist[1].get_tmrca(0, 1)
        data["pan"][r] = Ts
    with open("structure_simulated_results.pkl", "wb") as fout:
        pickle.dump(data, fout)
    return 


if not os.path.isfile("structure_simulated_results.pkl"):
    simulate()
with open("structure_simulated_results.pkl", "rb") as fin:
    data = pickle.load(fin)


fig, axes = plt.subplot_mosaic(
    "AABBCCCC;DDDDEEFF;.HHIIJJ.", figsize=(9, 7.5), layout="constrained")


demesdraw.tubes(g_struc, ax=axes["A"])
demesdraw.tubes(g_pan, ax=axes["B"])


axes["C"].plot(ts, str_rates, color="tab:blue", marker="x", label="structure model")
axes["C"].plot(ts, pan_rates, color="tab:orange", marker="o", label="panmictic model", markerfacecolor="none")
axes["C"].set_xlim(0,)
axes["C"].set_ylim(0,)
axes["C"].set_xlabel("$t$ b.p, generations")
axes["C"].set_ylabel("$\lambda(t)$")
axes["C"].legend(fontsize=8, framealpha=0)


y_struc = moments.Demes.LD(g_struc, ['S'], u=u, r=rs).H2(0)
y_pan = moments.Demes.LD(g_pan, ['S'], u=u, r=rs).H2(0)
axes["D"].plot(rs, y_struc, color="tab:blue", label="structure model")
axes["D"].plot(rs, y_pan, color="tab:orange", label="panmictic model")
axes["D"].set_xlabel("$r$")
axes["D"].set_ylabel("E[D^+]$")
axes["D"].set_xscale("log")
axes["D"].set_ylim(0,)
axes["D"].legend(fontsize=8, framealpha=0)


for label, ax in axes.items():
    ax.set_title(label, loc='left', fontsize='large')


axes["H"].scatter(data["struc"][rrs[0]][:, 0], data["struc"][rrs[0]][:, 1], 
    color="tab:blue", marker="o", s=4, linewidths=0)
axes["H"].scatter(data["pan"][rrs[0]][:, 0], data["pan"][rrs[0]][:, 1], 
    color="tab:orange", marker="o", s=4, linewidths=0)
axes["H"].set_ylim(0, 1e5)
axes["H"].set_xlim(0, 1e5)
axes["H"].set_ylabel("$T_y$")
axes["H"].set_xlabel("$T_x$")
axes["H"].set_title(r"$r=10^{-6}$")


axes["I"].scatter(data["struc"][rrs[1]][:, 0], data["struc"][rrs[1]][:, 1], 
    color="tab:blue", marker="o", s=4, linewidths=0)
axes["I"].scatter(data["pan"][rrs[1]][:, 0], data["pan"][rrs[1]][:, 1], 
    color="tab:orange", marker="o", s=4, linewidths=0)
axes["I"].set_ylim(0, 1e5)
axes["I"].set_xlim(0, 1e5)
axes["I"].set_ylabel("$T_y$")
axes["I"].set_xlabel("$T_x$")
axes["I"].set_title(r"$r=10^{-5}$")


axes["J"].scatter(data["struc"][rrs[2]][:, 0], data["struc"][rrs[2]][:, 1], 
    color="tab:blue", marker="o", s=4, linewidths=0)
axes["J"].scatter(data["pan"][rrs[2]][:, 0], data["pan"][rrs[2]][:, 1], 
    color="tab:orange", marker="o", s=4, linewidths=0)
axes["J"].set_ylim(0, 1e5)
axes["J"].set_xlim(0, 1e5)
axes["J"].set_ylabel("$T_y$")
axes["J"].set_xlabel("$T_x$")
axes["J"].set_title(r"$r=10^{-4}$")


plt.savefig("figure_theory_structure_idability.png", dpi=244)
