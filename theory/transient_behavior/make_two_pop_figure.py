
import numpy as np
import moments
import demes
import matplotlib.pyplot as plt 


def im_model(rs, t):
    b = demes.Builder()
    b.add_deme("A", epochs=[{"start_size": 10000, "end_time": t}])
    b.add_deme('popA', ancestors=["A"], epochs=[dict(end_time=0, start_size=10000)])
    b.add_deme('popB', ancestors=["A"], epochs=[dict(end_time=0, start_size=10000)])
    b.add_migration(demes=["popA", "popB"], rate=1e-4)
    g = b.resolve()
    y = moments.Demes.LD(g, ["popA", "popB"], u=u, r=rs)
    return y.H2(0), y.H()[0], y.H2(0, 1, phased=False), y.H()[1]


# Parameters
u = 1.4e-8
rs = np.logspace(-6, -2, 20)


fig, axes = plt.subplot_mosaic("ABC;DEF", figsize=(7, 5), layout="constrained")


# IM model
y0, H0, yij0, Hij0 = im_model(rs, 1e-5)
y1, H1, yij1, Hij1 = im_model(rs, 1e4)
y2, H2, yij2, Hij2 = im_model(rs, 2e5)

axes["A"].plot(rs, y0, color="tab:blue")
axes["A"].plot(rs, y1, color="tab:blue", linestyle="dotted")
axes["A"].plot(rs, y2, color="tab:blue")
axes["A"].plot(rs, yij0, color="tab:orange")
axes["A"].plot(rs, yij1, color="tab:orange", linestyle="dotted")
axes["A"].plot(rs, yij2, color="tab:orange")
axes["A"].set_xscale("log")
axes["A"].set_ylim(0, 2.6e-6)
axes["A"].set_ylabel(r"$E[D^+]$")
axes["A"].set_xlabel(r"$\rho$")

axes["C"].plot(rs, y0 / H0 ** 2, color="tab:blue")
axes["C"].plot(rs, y1 / H1 ** 2, color="tab:blue", linestyle="dotted")
axes["C"].plot(rs, y2 / H2 ** 2, color="tab:blue")
axes["C"].plot(rs, yij0 / Hij0 ** 2, color="tab:orange")
axes["C"].plot(rs, yij1 / Hij1 ** 2, color="tab:orange", linestyle="dotted")
axes["C"].plot(rs, yij2 / Hij2 ** 2, color="tab:orange")
axes["C"].set_xscale("log")
axes["C"].set_ylim(0.9,)
axes["C"].set_ylabel(r"$E[D^+] / E[H]^2$")
axes["C"].set_xlabel(r"$\rho$")

times = np.linspace(1e-5, 2e5, 20)
y0 = [im_model([1e-6], t) for t in times]
y1 = [im_model([1.8e-5], t) for t in times]
y2 = [im_model([3.6e-5], t) for t in times]
y3 = [im_model([1e-2], t) for t in times]
axes["B"].plot(times, [y[0] for y in y0], color="tab:blue")
axes["B"].plot(times, [y[0] for y in y1], color="tab:blue")
axes["B"].plot(times, [y[0] for y in y2], color="tab:blue")
axes["B"].plot(times, [y[0] for y in y3], color="tab:blue")
axes["B"].plot(times, [y[2] for y in y0], color="tab:orange")
axes["B"].plot(times, [y[2] for y in y1], color="tab:orange")
axes["B"].plot(times, [y[2] for y in y2], color="tab:orange")
axes["B"].plot(times, [y[2] for y in y3], color="tab:orange")
axes["B"].set_ylim(0, 2.6e-6)
axes["B"].set_xlim(0,)
axes["B"].set_ylabel("$E[D^+]$")
axes["B"].set_xlabel(r"$t$ since sep.")


# Pulse model
def pulse_model(rs, t):
    b = demes.Builder()
    if t <= 10000:
        b.add_deme("A", epochs=[{"start_size": 10000, "end_time": t}])
        b.add_deme('popA', ancestors=["A"], epochs=[dict(end_time=0, start_size=10000)])
        b.add_deme('popB', ancestors=["A"], epochs=[dict(end_time=0, start_size=1000)])
    else:
        b.add_deme("A", epochs=[{"start_size": 10000, "end_time": t}])
        b.add_deme('popA', ancestors=["A"], epochs=[dict(end_time=0, start_size=10000)])
        b.add_deme('popB', ancestors=["A"], epochs=[dict(end_time=0, start_size=1000)])
        b.add_pulse(sources=["popB"], dest="popA", proportions=[0.1], time=t - 10000)
    g = b.resolve()
    y = moments.Demes.LD(g, ["popA", "popB"], u=u, r=rs)
    return y.H2(0), y.H()[0], y.H2(0, 1, phased=False), y.H()[1]


y0, H0, yij0, Hij0 = pulse_model(rs, 9999)
y1, H1, yij1, Hij1 = pulse_model(rs, 10001)
y2, H2, yij2, Hij2 = pulse_model(rs, 2e4)

axes["D"].plot(rs, y0, color="tab:blue")
axes["D"].plot(rs, y1, color="tab:blue", linestyle="dotted")
axes["D"].plot(rs, y2, color="tab:blue")
axes["D"].plot(rs, yij0, color="tab:orange")
axes["D"].plot(rs, yij1, color="tab:orange", linestyle="dotted")
axes["D"].plot(rs, yij2, color="tab:orange")
axes["D"].set_xscale("log")
axes["D"].set_ylim(0,)
axes["D"].set_ylabel(r"$E[D^+]$")
axes["D"].set_xlabel(r"$\rho$")

axes["F"].plot(rs, y0 / H0 ** 2, color="tab:blue")
axes["F"].plot(rs, y1 / H1 ** 2, color="tab:blue", linestyle="dotted")
axes["F"].plot(rs, y2 / H2 ** 2, color="tab:blue")
axes["F"].plot(rs, yij0 / Hij0 ** 2, color="tab:orange")
axes["F"].plot(rs, yij1 / Hij1 ** 2, color="tab:orange", linestyle="dotted")
axes["F"].plot(rs, yij2 / Hij2 ** 2, color="tab:orange")
axes["F"].set_xscale("log")
axes["F"].set_ylim(0.9,)
axes["F"].set_ylabel(r"$E[D^+] / E[H]^2$")
axes["F"].set_xlabel(r"$\rho$")

times = np.linspace(1e-5, 2e4, 20)
y0 = [pulse_model([1e-6], t) for t in times]
y1 = [pulse_model([1.8e-5], t) for t in times]
y2 = [pulse_model([3.6e-5], t) for t in times]
y3 = [pulse_model([1e-2], t) for t in times]
axes["E"].plot(times, [y[0] for y in y0], color="tab:blue")
axes["E"].plot(times, [y[0] for y in y1], color="tab:blue")
axes["E"].plot(times, [y[0] for y in y2], color="tab:blue")
axes["E"].plot(times, [y[0] for y in y3], color="tab:blue")
axes["E"].plot(times, [y[2] for y in y0], color="tab:orange")
axes["E"].plot(times, [y[2] for y in y1], color="tab:orange")
axes["E"].plot(times, [y[2] for y in y2], color="tab:orange")
axes["E"].plot(times, [y[2] for y in y3], color="tab:orange")
axes["E"].set_ylim(0,)
axes["E"].set_xlim(0,)
axes["E"].set_ylabel("$E[D^+]$")
axes["E"].set_xlabel(r"$t$ since sep.")


for label, ax in axes.items():
    ax.set_title(label, loc='left', fontsize='large')


plt.savefig("figure_theory_two_pop.png", dpi=244)

