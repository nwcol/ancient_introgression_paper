
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




def equilibrium():
    """
    Plot D+ trajectories under a model of population expansion from 1e4 to 1e5.
    """
    u = 1.5e-8
    r_vals = np.logspace(-6, -2, 60)
    colors = palettes.Category10_10

    b = demes.Builder()
    b.add_deme('popA', epochs=[dict(end_time=0, start_size=10000)])
    g = b.resolve()
    y = moments.Demes.LD(g, ['popA'], u=u, r=r_vals)
    ys = y.LD()
    DP = y.H2(0)
    
    fig, axs = plt.subplots(1, 1, figsize=(4, 4), layout='constrained')
    ax = axs
    ax.plot(r_vals, DP, label='$D^+$', color=colors[0])
    ax.plot(r_vals, ys[:, 0], label='$D^2$', color=colors[1])
    ax.plot(r_vals, ys[:, 1], label='$Dz$', color=colors[2])
    ax.plot(r_vals, ys[:, 2], label='$\pi_2$', color=colors[3])
    ax.legend(framealpha=0)
    ax.set_xlabel('$r$')
    ax.set_xscale('log')
    ax.set_ylim(0,)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('figures/equilibrium_Dplus_decomp.png', dpi=244)
    plt.show()

    fig, axs = plt.subplots(1, 1, figsize=(4, 4), layout='constrained')
    ax = axs
    ax.plot(r_vals, DP, label='$D^+$', color=colors[0])
    ax.plot(r_vals, 4 * ys[:, 0], label='$4D^2$', color=colors[1])
    ax.plot(r_vals, 2 * ys[:, 1], label='$2Dz$', color=colors[2])
    ax.plot(r_vals, 4 * ys[:, 2], label='$4\pi_2$', color=colors[3])
    ax.legend(framealpha=0)
    ax.set_xlabel('$r$')
    ax.set_xscale('log')
    ax.set_ylim(0,)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig('figures/equilibrium_Dplus_decomp_scaled.png', dpi=244)
    plt.show()

    return 


equilibrium()
