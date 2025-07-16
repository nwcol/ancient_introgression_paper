
import numpy as np
import moments
import demes
import matplotlib.pyplot as plt 


def expansion_model(rs, t):
    b = demes.Builder()
    b.add_deme('popA', epochs=[
        dict(end_time=t, start_size=10000), dict(end_time=0, start_size=20000)])
    g = b.resolve()
    y = moments.Demes.LD(g, ['popA'], u=u, r=rs)
    return y.H2(0), y.H()[0]


# Parameters
u = 1.4e-8
Ne0 = 10000
rs = np.logspace(-6, -2, 20)


fig, axes = plt.subplot_mosaic("ABC;DEF", figsize=(7, 5), layout="constrained")


# Twofold expansion model
y0, H0 = expansion_model(rs, 1e-5)
y1, H1 = expansion_model(rs, 2e4)
y2, H2 = expansion_model(rs, 1e6)

axes["A"].plot(rs, y0, color="black")
axes["A"].plot(rs, y1, color="black", linestyle="dotted")
axes["A"].plot(rs, y2, color="black")
axes["A"].set_xscale("log")
axes["A"].set_ylim(0, 2.6e-6)
axes["A"].set_ylabel(r"$E[D^+]$")
axes["A"].set_xlabel(r"$\rho$")

axes["C"].plot(rs, y0 / H0 ** 2, color="black")
axes["C"].plot(rs, y1 / H1 ** 2, color="black", linestyle="dotted")
axes["C"].plot(rs, y2 / H2 ** 2, color="black")
axes["C"].set_xscale("log")
axes["C"].set_ylim(0.9, )
axes["C"].set_ylabel(r"$E[D^+] / E[H]^2$")
axes["C"].set_xlabel(r"$\rho$")

times = np.linspace(1e-5, 2e5, 20)
y0 = np.array([expansion_model([1e-6], t)[0] for t in times])
y1 = np.array([expansion_model([1.8e-5], t)[0] for t in times])
y2 = np.array([expansion_model([3.6e-5], t)[0] for t in times])
y3 = np.array([expansion_model([1e-2], t)[0] for t in times])
axes["B"].plot(times, y0, color="black")
axes["B"].plot(times, y1, color="black")
axes["B"].plot(times, y2, color="black")
axes["B"].plot(times, y3, color="black")
axes["B"].set_ylim(0, 2.6e-6)
axes["B"].set_xlim(0, )
axes["B"].set_ylabel("$E[D^+]$")
axes["B"].set_xlabel(r"$t$ since exp.")


# Twofold contraction model
def contraction_model(rs, t):
    b = demes.Builder()
    b.add_deme('popA', epochs=[
        dict(end_time=t, start_size=10000), dict(end_time=0, start_size=5000)])
    g = b.resolve()
    y = moments.Demes.LD(g, ['popA'], u=u, r=rs)
    return y.H2(0), y.H()[0]

y0, H0 = contraction_model(rs, 1e-5)
y1, H1 = contraction_model(rs, 1e4)
y2, H2 = contraction_model(rs, 1e6)

axes["D"].plot(rs, y0, color="black")
axes["D"].plot(rs, y1, color="black", linestyle="dotted")
axes["D"].plot(rs, y2, color="black")
axes["D"].set_xscale("log")
axes["D"].set_ylim(0, 6.5e-7)
axes["D"].set_ylabel(r"$E[D^+]$")
axes["D"].set_xlabel(r"$\rho$")

axes["F"].plot(rs, y0 / H0 ** 2, color="black")
axes["F"].plot(rs, y1 / H1 ** 2, color="black", linestyle="dotted")
axes["F"].plot(rs, y2 / H2 ** 2, color="black")
axes["F"].set_xscale("log")
axes["F"].set_ylim(0.9,)
axes["F"].set_ylabel(r"$E[D^+] / E[H]^2$")
axes["F"].set_xlabel(r"$\rho$")

times = np.linspace(1, 5e4, 20)
y0 = np.array([contraction_model([1e-6], t)[0] for t in times])
y1 = np.array([contraction_model([3.6e-5], t)[0] for t in times])
y2 = np.array([contraction_model([7.2e-5], t)[0] for t in times])
y3 = np.array([contraction_model([1e-2], t)[0] for t in times])
axes["E"].plot(times, y0, color="black")
axes["E"].plot(times, y1, color="black")
axes["E"].plot(times, y2, color="black")
axes["E"].plot(times, y3, color="black")
axes["E"].set_ylim(0, 6.5e-7)
axes["E"].set_xlim(0,)
axes["E"].set_ylabel("$E[D^+]$")
axes["E"].set_xlabel(r"$t$ since contr.")


for label, ax in axes.items():
    ax.set_title(label, loc='left', fontsize='large')


plt.savefig("figure_theory_one_pop.png", dpi=244)

