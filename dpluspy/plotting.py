
from bokeh import palettes
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy import stats

from . import parsing, inference, utils


mpl.rcParams["xtick.labelsize"] = 8
mpl.rcParams["ytick.labelsize"] = 8
mpl.rcParams["font.size"] = 10
mpl.rcParams["axes.titlesize"] = 10
mpl.rcParams["font.style"] = "normal"
mpl.rcParams["font.family"] = "sans-serif"
mpl.rcParams["savefig.bbox"] = "tight"



def plot_d_plus_curves(
    stats,
    stats_to_plot=[],
    rows=None,
    cols=None,
    ax_size=2,
    dpi=244,
    bins=None,
    rs=None,
    cM=False,
    out_file=None,
    show=True,
):
    """
    
    """
    statistics = stats.names()
    labels = utils.get_latex_names(stats.pop_ids)
    if len(stats_to_plot) == 0:
        stats_to_plot = stats.names()[0]
    num_stats = len(stats_to_plot)
    if not cols:
        cols = min(5, num_stats)
    if not rows:
        rows = -(num_stats // -cols)
    figsize = (cols * ax_size, rows * ax_size)
    bins = np.asarray(bins)
    mids = (bins[1:] + bins[:-1]) / 2
    x_label = "$r$"
    if cM == True:
        mids = utils.map_function(mids)
        x_label = "cM"

    fig, axs = plt.subplots(rows, cols, figsize=figsize, layout="constrained")
    axs = axs.flat
    for ax in axs[num_stats:]:
        ax.remove()

    for i, stat in enumerate(stats_to_plot):
        ax = axs[i]
        k = statistics[0].index(stat)
        y = [stats[j][k] for j in range(len(mids))]
        ax.plot(mids, y)
        ax.set_xscale("log")
        if i >= (cols * rows - rows):
            ax.set_xlabel(x_label)
        if i % cols == 0:
            ax.set_ylabel("$D^+$")
        label = labels[k]
        ax.set_title(label, y=0.85)

    if out_file:
        plt.savefig(out_file, dpi=dpi)
    if show:
        fig.show()

    return fig


def plot_d_plus_curves_comparison(
    model,
    means,
    varcovs,
    stats_to_plot=[],
    fill=True,
    rows=None,
    cols=None,
    ax_size=2,
    dpi=244,
    bins=None,
    rs=None,
    cM=False,
    out_file=None,
    show=True
):
    """
    
    """
    statistics = model.names()
    labels = utils.get_latex_names(model.pop_ids)
    if len(stats_to_plot) == 0:
        stats_to_plot = model.names()[0]
    num_stats = len(stats_to_plot)
    if not cols:
        cols = min(5, num_stats)
    if not rows:
        rows = -(num_stats // -cols)
    figsize = (cols * ax_size, rows * ax_size)

    bins = np.asarray(bins)
    mids = (bins[1:] + bins[:-1]) / 2
    x_label = "$r$"
    if cM == True:
        mids = utils.map_function(mids)
        x_label = "cM"

    fig, axs = plt.subplots(rows, cols, figsize=figsize, layout="constrained")
    axs = axs.flat
    for ax in axs[num_stats:]:
        ax.remove()

    color = palettes.Category10_10[0]

    for i, stat in enumerate(stats_to_plot):
        ax = axs[i]
        k = statistics[0].index(stat)
        y_err = np.array(
            [varcovs[j][k, k] ** 0.5 * 1.96 for j in range(len(mids))]
        )
        y_data = np.array([means[j][k] for j in range(len(mids))])
        ax.fill_between(mids, y_data - y_err, y_data + y_err, alpha=0.30)
        ax.plot(mids, y_data, linestyle="dotted", color=color)

        y_exp = [model[j][k] for j in range(len(mids))]
        ax.plot(mids, y_exp, color=color)
        ax.set_xscale("log")
        if i >= (cols * rows - rows):
            ax.set_xlabel(x_label)
        if i % cols == 0:
            ax.set_ylabel("$D^+$")
        label = labels[k]
        ax.set_title(label, y=0.85)

    if out_file:
        plt.savefig(out_file, dpi=dpi)
    if show:
        fig.show()

    return fig


def plot_empirical_d_plus_curves(
    means,
    varcovs,
    pop_ids,
    stats_to_plot=[],
    pops_to_plot=[],
    fill=True,
    rows=None,
    cols=None,
    ax_size=2,
    dpi=244,
    bins=None,
    rs=None,
    cM=False,
    out_file=None,
    show=True,
    labels=None
):
    """
    
    """
    statistics = utils.stat_names(pop_ids)
    stat_labels = utils.get_latex_names(pop_ids)
    if len(stats_to_plot) == 0 and len(pops_to_plot) == 0:
        stats_to_plot = statistics[0]
    elif len(stats_to_plot) == 0 and len(pops_to_plot) > 0:
        sorted_pops = [pop for pop in pop_ids if pop in pops_to_plot]
        stats_to_plot = utils.stat_names(sorted_pops)[0]
    num_stats = len(stats_to_plot)
    if not cols:
        cols = min(5, num_stats)
    if not rows:
        rows = -(num_stats // -cols)
    figsize = (cols * ax_size, rows * ax_size)

    bins = np.asarray(bins)
    mids = (bins[1:] + bins[:-1]) / 2
    x_label = "$r$"
    if cM == True:
        mids = utils.map_function(mids)
        x_label = "cM"

    fig, axs = plt.subplots(rows, cols, figsize=figsize, layout="constrained")
    if rows > 1:
        axs = axs.flat
    elif cols == 1:
        axs = [axs]
    for ax in axs[num_stats:]:
        ax.remove()

    if not isinstance(means[0], list):
        means = [means]
        varcovs = [varcovs]

    for k, (_means, _varcovs) in enumerate(zip(means, varcovs)):
        label = labels[k] if labels is not None else None
        color = palettes.Category10_10[k]
        for i, stat in enumerate(stats_to_plot):
            ax = axs[i]
            l = statistics[0].index(stat)
            y_err = np.array(
                [_varcovs[j][l, l] ** 0.5 * 1.96 for j in range(len(mids))]
            )
            y_data = np.array([_means[j][l] for j in range(len(mids))])
            if fill:
                ax.fill_between(
                    mids, y_data - y_err, y_data + y_err, alpha=0.30
                )
                ax.plot(mids, y_data, linestyle="dotted", color=color, label=label)
            else:
                ax.errorbar(
                    mids, 
                    y_data, 
                    yerr=y_err, 
                    capsize=0, 
                    markersize=4,
                    elinewidth=1,
                    fmt="o",
                    mfc="none",
                    mec=color,
                    markeredgewidth=1,
                    ecolor=color,
                    label=label
                )
            if label is not None:
                ax.legend()
            ax.set_xscale("log")
            if i >= (cols * rows - rows):
                ax.set_xlabel(x_label)
            if i % cols == 0:
                ax.set_ylabel("$D^+$")
            stat_label = stat_labels[i]
            ax.set_title(stat_label, y=0.85)

    if out_file:
        plt.savefig(out_file, dpi=dpi)
    if show:
        fig.show()

    return fig
