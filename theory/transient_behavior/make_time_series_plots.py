
from bokeh import palettes
import demes
import matplotlib.pyplot as plt
import matplotlib as mpl
import moments
import numpy as np


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = r"\usepackage{amsmath}\usepackage{amssymb}"
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.serif'] = "Computer Modern"
mpl.rcParams['savefig.bbox'] = "tight"


u = 1.5e-8


def expansion():
    """
    Plot D+ trajectories under a model of population expansion from 1e4 to 1e5.
    """
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    ax0, ax1, ax2 = axs

    r_vals = np.logspace(-6, -2, 20)
    t_vals = [1, 30000, 1e6]
    labels = ['$t=0$', '$t=3 \cdot 10^4$', '$t=\infty$']
    styles = ['solid', 'dotted', 'dashed', 'dashdot']

    ys = []
    for t in t_vals:
        b = demes.Builder()
        b.add_deme('popA', epochs=[dict(end_time=t, start_size=10000), 
                                dict(end_time=0, start_size=20000)])
        g = b.resolve()
        y = moments.Demes.LD(g, ['popA'], u=u, r=r_vals).H2(0)
        ys.append(y)
    for i in range(3):
        ax0.plot(r_vals, ys[i], label=labels[i], color='black', linestyle=styles[i])
    ax0.legend(framealpha=0)
    ax0.set_xlabel('$r$')
    ax0.set_ylabel('$E[D^+]$')
    ax0.set_xscale('log')
    ax0.set_ylim(0,)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)

    
    for i in range(3):
        ax1.plot(r_vals, ys[i]/ys[i][-1], label=labels[i], color='black', linestyle=styles[i])
    ax1.legend(framealpha=0)
    ax1.set_xlabel('$r$')
    ax1.set_ylabel('E[$D^+/H^2]$')
    ax1.set_xscale('log')
    ax1.set_ylim(0.9,)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)


    r_vals = np.array([1e-6, 1.8e-5, 3.6e-5, 1e-2])
    t_vals = np.linspace(1, 200000, 50)
    ys = []
    for t in t_vals:
        b = demes.Builder()
        b.add_deme('popA', epochs=[dict(end_time=t, start_size=10000), 
                                dict(end_time=0, start_size=20000)])
        g = b.resolve()
        y = moments.Demes.LD(g, ['popA'], u=u, r=r_vals).H2(0)
        ys.append(y)
    ys = np.array(ys)
    ax2.plot(t_vals, ys[:, 0], label='$r=10^{-6}$', color='black')
    ax2.plot(t_vals, ys[:, 1], label='$r=4 \cdot 10^{-6}$', color='black', linestyle='dashed')
    ax2.plot(t_vals, ys[:, 2], label='$r=4 \cdot 10^{-5}$', color='black', linestyle='dotted')
    ax2.plot(t_vals, ys[:, 3], label='$r=10^{-2}$', color='black', linestyle='dashdot')
    ax2.legend(framealpha=0)
    ax2.set_xlabel('time since two-fold expansion (generations)')
    ax2.set_ylabel('$E[D^+]$')
    ax2.set_xlim(0,)
    ax2.set_ylim(0,)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # ax.invert_xaxis()


    plt.savefig('figure_expansion_trajectories.png', dpi=244)
    plt.show()

    return 


def contraction():
    """
    Plot D+ trajectories under a model of contraction from 1e4 to 1e3
    """
    fig, axs = plt.subplots(2, 2, figsize=(7, 7), layout='constrained')
    ax0, ax1, ax2, ax3 = axs.flat

    r_vals = np.logspace(-6, -2, 20)
    t_vals = [1, 15000, 1e6]
    labels = ['$t=0$', '$t=1.5 \cdot 10^4$', '$t=\infty$']
    styles = ['solid', 'dotted', 'dashed', 'dashdot']

    ys = []
    hs = []
    for t in t_vals:
        b = demes.Builder()
        b.add_deme('popA', epochs=[dict(end_time=t, start_size=10000), 
                                dict(end_time=0, start_size=5000)])
        g = b.resolve()
        model = moments.Demes.LD(g, ['popA'], u=u, r=r_vals)
        ys.append(model.H2(0))
        hs.append(model.H(pops=[0]))
    ys = np.array(ys)
    hs = np.array(hs)

    # subfig A
    for i in range(3):
        ax0.plot(r_vals, ys[i], label=labels[i], color='black', linestyle=styles[i])
    ax0.legend(framealpha=0)
    ax0.set_xlabel('$r$')
    ax0.set_ylabel('$E[D^+]$')
    ax0.set_xscale('log')
    ax0.set_ylim(0,7.6e-7)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)

    # subfig C
    for i in range(3):
        ax2.plot(r_vals, ys[i]/hs[i]**2, label=labels[i], color='black', linestyle=styles[i])
    ax2.legend(framealpha=0)
    ax2.set_xlabel('$r$')
    ax2.set_ylabel('E[$D^+]/E[\pi]^2$')
    ax2.set_xscale('log')
    ax2.set_ylim(0.8,3.1)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    ax2.grid(axis='y', alpha=0.2)

    r_vals = np.array([1e-6, 3.6e-5, 7.2e-5, 1e-2])
    labels = ['$r=10^{-6}$', '$r=3.6 \cdot 10^{-5}$', '$r=7.2 \cdot 10^{-5}$', '$r=10^{-2}$']
    t_vals = np.linspace(1, 50000, 40)
    ys = []
    hs = []
    for t in t_vals:
        b = demes.Builder()
        b.add_deme('popA', epochs=[dict(end_time=t, start_size=10000), 
                                dict(end_time=0, start_size=5000)])
        g = b.resolve()
        model = moments.Demes.LD(g, ['popA'], u=u, r=r_vals)
        ys.append(model.H2(0))
        hs.append(model.H(pops=[0]))
    ys = np.array(ys)
    hs = np.array(hs)

    # subfig B
    for i in range(4):
        ax1.plot(t_vals, ys[:, i], label=labels[i], color='black', linestyle=styles[i])
    ax1.legend(framealpha=0)
    ax1.set_xlabel('time since two-fold contraction (generations)')
    ax1.set_ylabel('$E[D^+]$')
    ax1.set_xlim(0,)
    ax1.set_ylim(0,7.6e-7)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)

    # subfig D
    for i in range(4):
        ax3.plot(t_vals, ys[:, i]/hs[:, 0]**2, label=labels[i], color='black', linestyle=styles[i])
    ax3.legend(framealpha=0)
    ax3.set_xlabel('time since two-fold contraction (generations)')
    ax3.set_ylabel('E[$D^+]/E[\pi]^2$')
    ax3.set_xlim(0,)
    ax3.set_ylim(0.8,3.1)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)

    plt.savefig('figure_contraction_trajectories.png', dpi=244)
    plt.show()

    return


def bottleneck():
    """
    Plot D+ trajectories under a model of bottleneck. At t=0 the population
    transitions from Ne=1e4 to 1e3, then at t=2000 recovers to 1e4. We plot
    snapshots following the recovery also.
    """
    fig, axs = plt.subplots(1, 3, figsize=(12, 4), layout='constrained')
    ax0, ax1, ax2 = axs

    r_vals = np.logspace(-6, -2, 20)
    t_vals = np.array([1, 2000, 20000, 1e6])
    labels = ['$t=0$', '$t=2 \cdot 10^3$', '$t=2 \cdot 10^4$', '$t=\infty$']
    styles = ['solid', 'dotted', 'dashed', 'dashdot']

    ys = []
    for t in t_vals:
        b = demes.Builder()
        b.add_deme('popA', epochs=[dict(end_time=t + 2000, start_size=10000),
                                   dict(end_time=t, start_size=1000), 
                                   dict(end_time=0, start_size=10000)])
        g = b.resolve()
        y = moments.Demes.LD(g, ['popA'], u=u, r=r_vals).H2(0)
        ys.append(y)
    for i in range(4):
        ax0.plot(r_vals, ys[i], label=labels[i], color='black', linestyle=styles[i])
    ax0.legend(framealpha=0)
    ax0.set_xlabel('$r$')
    ax0.set_ylabel('$E[D^+]$')
    ax0.set_xscale('log')
    ax0.set_ylim(0,)
    ax0.spines['top'].set_visible(False)
    ax0.spines['right'].set_visible(False)

    
    for i in range(4):
        ax1.plot(r_vals, ys[i]/ys[i][-1], label=labels[i], color='black', linestyle=styles[i])
    ax1.legend(framealpha=0)
    ax1.set_xlabel('$r$')
    ax1.set_ylabel('E[$D^+/H^2]$')
    ax1.set_xscale('log')
    ax1.set_ylim(0.9,)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)


    # plot with xax t, yax D+, key r
    r_vals = np.array([1e-6, 3.6e-5, 3.6e-4, 1e-2])
    labels = ['$r=10^{-6}$', '$r=3.6 \cdot 10^{-6}$', '$r=3.6 \cdot 10^{-4}$', '$r=10^{-2}$']
    t_vals = np.linspace(1, 100000, 50)
    ys = []
    for t in t_vals:
        b = demes.Builder()
        b.add_deme('popA', epochs=[dict(end_time=t + 2000, start_size=10000),
                                   dict(end_time=t, start_size=1000), 
                                   dict(end_time=0, start_size=10000)])
        g = b.resolve()
        y = moments.Demes.LD(g, ['popA'], u=u, r=r_vals).H2(0)
        ys.append(y)
    ys = np.array(ys)
    for i in range(4):
        ax2.plot(t_vals, ys[:, i], label=labels[i], color='black', linestyle=styles[i])
    ax2.legend(framealpha=0)
    ax2.set_xlabel('time since bottleneck recovery (generations)')
    ax2.set_ylabel('$E[D^+]$')
    ax2.set_xlim(0,)
    ax2.set_ylim(0,)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)
    # ax.invert_xaxis()


    plt.savefig('figure_bottleneck_trajectories.png', dpi=244)
    plt.show()

    return


def pulsed_introgression_gamma05():
    """
    Plot D+ trajectories under a model of pulsed admixture.
    """
    u = 1.5e-8

    def model_func(_t, _rs):
        b = demes.Builder()
        if _t < 0:
            b.add_deme('anc', epochs=[dict(end_time=20000 + _t, start_size=10000)])
            b.add_deme('popA', ancestors=['anc'], 
                epochs=[dict(end_time=0, start_size=10000)])
            b.add_deme('popB', ancestors=['anc'],
                epochs=[dict(end_time=0, start_size=10000)])
        else:
            b.add_deme('anc', epochs=[dict(end_time=20000 + _t, start_size=10000)])
            b.add_deme('popA', ancestors=['anc'], 
                epochs=[dict(end_time=0, start_size=10000)])
            b.add_deme('popB', ancestors=['anc'],
                epochs=[dict(end_time=0, start_size=10000)])
            b.add_pulse(sources=['popB'], proportions=[0.10], dest='popA', time=_t)
        g = b.resolve()
        model = moments.Demes.LD(g, ['popA', 'popB'], u=u, r=_rs)
        dA = model.H2(0)
        hA, hAB = model.H(pops=[0, 1])[:2]
        dAB = model.H2(0, 1, phased=False)
        return dA, hA, dAB, hAB

    fig, axs = plt.subplots(4, 2, figsize=(6, 12), layout='constrained')
    ax0, ax1, ax2, ax3, ax4, ax5, ax6, ax7 = axs.flat

    styles = ['solid', 'dotted', 'dashed', 'dashdot']   

    r_vals = np.array([1e-6, 3.6e-5, 7.2e-5, 1e-2])
    labels = ['$r=10^{-6}$', '$r=3.6 \cdot 10^{-5}$', '$r=7.2 \cdot 10^{-5}$', '$r=10^{-2}$']
    t_vals = np.linspace(-1000, 2000, 30)
    dAs, hAs, dABs, hABs = [], [], [], []
    for t in t_vals:
        dA, hA, dAB, hAB = model_func(t, r_vals)
        dAs.append(dA)
        hAs.append(hA)
        dABs.append(dAB)
        hABs.append(hAB)
    dAs, hAs, dABs, hABs = [np.array(x) for x in [dAs, hAs, dABs, hABs]]

    # Panel 1: D+_A time trajectory
    for i in range(4):
        ax1.plot(t_vals, dAs[:, i], label=labels[i], color='black', 
            linestyle=styles[i])
    ax1.set_xlabel('time since admixture (generations)')
    ax1.set_ylabel('$E[D^+_{A}]$')
    ax1.set_ylim(0,)

    # Panel 3: D+_A/H_A time trajectory
    for i in range(4):
        ax3.plot(t_vals, dAs[:, i] / hAs ** 2, label=labels[i], 
            color='black', linestyle=styles[i])
    ax3.set_xlabel('time since two-fold contraction (generations)')
    ax3.set_ylabel('E[$D^+_{A}]/E[H_{A}]^2$')
    ax3.set_ylim(0.8,)

    # Panel 5: D+_AB time trajectory
    for i in range(4):
        ax5.plot(t_vals, dABs[:, i], label=labels[i], color='black', 
            linestyle=styles[i])
    ax5.set_xlabel('time since admixture (generations)')
    ax5.set_ylabel('$E[D^+_{AB}]$')
    ax5.set_ylim(0,)

    # Panel 7: D+_AB/H_AB time trajectory
    for i in range(4):
        ax7.plot(t_vals, dABs[:, i] / hABs ** 2, label=labels[i], 
            color='black', linestyle=styles[i])
    ax7.set_xlabel('time since two-fold contraction (generations)')
    ax7.set_ylabel('E[$D^+_{AB}]/E[H_{AB}]^2$')
    ax7.set_ylim(0.8,)

    for ax in axs.flat:
        ax.legend(framealpha=0)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.savefig('figure_introgression_gamma10_trajectory.png', dpi=244)
    plt.show()

    return



pulsed_introgression_gamma05()


