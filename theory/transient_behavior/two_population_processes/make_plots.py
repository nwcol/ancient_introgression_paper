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



def pure_isolation():
    """
    
    """
    def model_func(_t, _rs):

        b = demes.Builder()
        b.add_deme('anc', epochs=[dict(end_time=_t, start_size=10000)])
        b.add_deme('popA', ancestors=['anc'], 
            epochs=[dict(end_time=0, start_size=5000)])
        b.add_deme('popB', ancestors=['anc'],
            epochs=[dict(end_time=0, start_size=5000)])
        g = b.resolve()
        model = moments.Demes.LD(g, ['popA', 'popB'], u=u, r=_rs)
        D_AB = model.H2(0, 1, phased=False)
        H_AB = model.H(pops=[0, 1])[1]
        return D_AB, H_AB

    u = 1.5e-8 

    fig, axs = plt.subplots(2, 2, figsize=(6, 6), layout='constrained')
    ax0, ax1, ax2, ax3 = axs.flat

    r_vals = np.logspace(-6, -2, 20)
    t_vals = [1, 20000, 40000]
    labels = ['$t=0$', '$t=3 \cdot 10^4$', '$t=\infty$']
    labels = ['$t=0$', '$t=3 \cdot 10^4$', '$t=\infty$']
    styles = ['solid', 'dotted', 'dashed', 'dashdot']
    ys = []
    hs = []
    for t in t_vals:
        y, h = model_func(t, r_vals)
        ys.append(y)
        hs.append(h)
    ys = np.array(ys)
    hs = np.array(hs)

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
        ax2.plot(r_vals, ys[i] / hs[i] ** 2, label=labels[i], color='black', linestyle=styles[i])
    ax2.legend(framealpha=0)
    ax2.set_xlabel('$r$')
    ax2.set_ylabel('$E[D^+]/E[\pi_{ij}]^2$')
    ax2.set_xscale('log')
    ax2.set_ylim(0,)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    r_vals = np.array([1e-6, 1.8e-5, 3.6e-5, 1e-2])
    t_vals = np.linspace(1, 40000, 50)
    time_ys = []
    time_hs = []
    for t in t_vals:
        y, h = model_func(t, r_vals)
        time_ys.append(y)
        time_hs.append(h)
    time_ys = np.array(time_ys)
    time_hs = np.array(time_hs)

    labels = ['$r=10^{-6}$', '$r=1.8 \cdot 10^{-5}$', '$r=3.6 \cdot 10^{-5}$', '$r=10^{-2}$']
    for i in range(4):
        ax1.plot(t_vals, time_ys[:, i], label=labels[i], color='black', linestyle=styles[i])
    ax1.legend(framealpha=0)
    ax1.set_xlabel('time since isolation (generations)')
    ax1.set_ylabel('$E[D^+]$')
    ax1.set_xlim(0,)
    ax1.set_ylim(0,)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.invert_xaxis()

    for i in range(4):
        ax3.plot(t_vals, time_ys[:, i] / time_hs[i] ** 2, label=labels[i], color='black', linestyle=styles[i])
    ax3.legend(framealpha=0)
    ax3.set_xlabel('time since isolation (generations)')
    ax3.set_ylabel('$E[D^+]/E[\pi_{ij}]^2$')
    ax3.set_xlim(0,)
    ax3.set_ylim(0,)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.invert_xaxis()
    
    fig.suptitle("Pure isolation")

    plt.savefig('figures/figure_pure_isolation.png', dpi=244)
    plt.show()

    return


def isolation_migration(m):
    """
    
    """
    def model_func(_t, _rs):

        b = demes.Builder()
        b.add_deme('anc', epochs=[dict(end_time=_t, start_size=10000)])
        b.add_deme('popA', ancestors=['anc'], 
            epochs=[dict(end_time=0, start_size=5000)])
        b.add_deme('popB', ancestors=['anc'],
            epochs=[dict(end_time=0, start_size=5000)])
        b.add_migration(demes=['popA', 'popB'], rate=m)
        g = b.resolve()
        model = moments.Demes.LD(g, ['popA', 'popB'], u=u, r=_rs)
        D_AB = model.H2(0, 1, phased=False)
        H_AB = model.H(pops=[0, 1])[1]
        return D_AB, H_AB

    u = 1.5e-8 

    fig, axs = plt.subplots(2, 2, figsize=(6, 6), layout='constrained')
    ax0, ax1, ax2, ax3 = axs.flat

    r_vals = np.logspace(-6, -2, 20)
    t_vals = [1, 20000, 40000]
    labels = ['$t=0$', '$t=2 \cdot 10^4$', '$t=4 \cdot 10^4$']
    styles = ['solid', 'dotted', 'dashed', 'dashdot']
    ys = []
    hs = []
    for t in t_vals:
        y, h = model_func(t, r_vals)
        ys.append(y)
        hs.append(h)
    ys = np.array(ys)
    hs = np.array(hs)

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
        ax2.plot(r_vals, ys[i] / hs[i] ** 2, label=labels[i], color='black', linestyle=styles[i])
    ax2.legend(framealpha=0)
    ax2.set_xlabel('$r$')
    ax2.set_ylabel('$E[D^+]/E[\pi_{ij}]^2$')
    ax2.set_xscale('log')
    ax2.set_ylim(0,)
    ax2.spines['top'].set_visible(False)
    ax2.spines['right'].set_visible(False)

    r_vals = np.array([1e-6, 1.8e-5, 3.6e-5, 1e-2])
    t_vals = np.linspace(1, 100000, 50)
    time_ys = []
    time_hs = []
    for t in t_vals:
        y, h = model_func(t, r_vals)
        time_ys.append(y)
        time_hs.append(h)
    time_ys = np.array(time_ys)
    time_hs = np.array(time_hs)

    labels = ['$r=10^{-6}$', '$r=1.8 \cdot 10^{-5}$', '$r=3.6 \cdot 10^{-5}$', '$r=10^{-2}$']
    for i in range(4):
        ax1.plot(t_vals, time_ys[:, i], label=labels[i], color='black', linestyle=styles[i])
    ax1.legend(framealpha=0)
    ax1.set_xlabel('time since isolation (generations)')
    ax1.set_ylabel('$E[D^+]$')
    ax1.set_xlim(0,)
    ax1.set_ylim(0,)
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.invert_xaxis()

    for i in range(4):
        ax3.plot(t_vals, time_ys[:, i] / time_hs[i] ** 2, label=labels[i], color='black', linestyle=styles[i])
    ax3.legend(framealpha=0)
    ax3.set_xlabel('time since isolation (generations)')
    ax3.set_ylabel('$E[D^+]/E[\pi_{ij}]^2$')
    ax3.set_xlim(0,)
    ax3.set_ylim(0,)
    ax3.spines['top'].set_visible(False)
    ax3.spines['right'].set_visible(False)
    ax3.invert_xaxis()

    fig.suptitle(f"Isolation with migration, $m={m}$")

    plt.savefig(f'figures/figure_isolation_migration_m{m}.png', dpi=244)
    plt.show()

    return


pure_isolation()
isolation_migration(1e-6)
isolation_migration(1e-5)
isolation_migration(1e-4)
isolation_migration(1e-3)
