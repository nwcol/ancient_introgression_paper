
import argparse
from bokeh.palettes import Category10_10, TolRainbow23
import demes
import matplotlib as mpl
import matplotlib.pyplot as plt
import moments
import numpy as np
import pickle

from h2py import parsing


def get_args():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-g', '--graph_file', type=str, required=True
    )
    parser.add_argument(
        '-d', '--data_file', type=str, default=None
    )
    parser.add_argument(
        '--pops', nargs='*', default=None, type=str,
    )
    parser.add_argument(
        '--cross_pop', action="store_true" 
    )
    parser.add_argument(
        '-o', '--out_file', required=True
    )
    parser.add_argument(
        '-u', '--u', type=float, required=True
    )
    parser.add_argument(
        '-T0', '--T0', type=float, default=None
    )
    parser.add_argument(
        '-n', '--n', type=int, default=None
    )
    parser.add_argument(
        '-log_time', '--log_time', action="store_true"
    )
    parser.add_argument(
        '-log_y', '--log_y', action="store_true"
    )
    parser.add_argument(
        '-title', '--title', type=str, default=None
    )

    return parser.parse_args()


def get_h2(g, t, u, rs, cross_pop=False):
    """
    Get expected D+ for every deme existing at time `t`.
    """
    deme_map = {deme.name: deme for deme in g.demes}
    sample_map = {}
    sample_pops = []
    pops = []
    ts = []
    b = demes.Builder.fromdict(g.asdict())
    for deme in g.demes:
        pop = deme.name
        if deme.epochs[-1].end_time <= t:
            if deme.epochs[0].start_time > t:
                t0 = t
            elif deme.epochs[0].start_time == t:
                t0 = t - 0.01
            else:
                continue
            sample_name = f"{pop}_SAMPLE"
            Ne = deme_map[pop].size_at(t0)
            b.add_deme(
                sample_name,
                ancestors=[pop],
                start_time=t0,
                epochs=[dict(end_time=t0 - 0.01, start_size=Ne)]
            )
            sample_pops.append(sample_name)
            pops.append(pop)
            ts.append(t0 - 0.01)
    g = b.resolve()
    name_map = {sample: pop for sample, pop in zip(sample_pops, pops)}
    print(t, ts, pops)
    y = moments.Demes.LD(g, sample_pops, sample_times=ts, r=rs, u=u)
    h2s = {}
    for i, pop in enumerate(sample_pops):
        h2s[name_map[pop]] = y.H2(i)
    if cross_pop:
        for i, pop_i in enumerate(sample_pops):
            for j, pop_j in enumerate(sample_pops):
                if i > j:
                    if name_map[pop_i] in deme_map[name_map[pop_j]].ancestors:
                        pass 
                    elif name_map[pop_j] in deme_map[name_map[pop_i]].ancestors:
                        pass
                    else:
                        h2s[(name_map[pop_i], name_map[pop_j])] = y.H2(i, j, phased=False)

    return h2s


def load_data(g, rs, data_file):

    demes = [deme.name for deme in g.demes]
    with open(data_file, "rb+") as f:
        data_dic = pickle.load(f)
    data = data_dic[next(iter(data_dic))]
    pops = data["pops"]
    pops_to_keep = [pop for pop in pops if pop in demes]
    deme_data = {}
    _data = parsing.subset_statistics(data, to_pops=pops_to_keep)
    bins = _data["bins"] 
    midpoints = bins[1:] + (bins[1:] - bins[:-1]) / 2
    labels = []
    for i in range(len(pops_to_keep)):
        for j in range(i, len(pops_to_keep)):
            if i == j:
                labels.append(pops_to_keep[i])
            else:
                labels.append((pops_to_keep[i], pops_to_keep[j]))
    for k, label in enumerate(labels):
        means = _data["means"][:-1, k]
        approx_h2 = np.interp(rs, midpoints, means, left=means[0])
        deme_data[label] = approx_h2
    print(deme_data.keys())
    return deme_data


def main():

    args = get_args()
    g = demes.load(args.graph_file)
    if args.T0 is not None:
        t_max = args.T0
    else:
        # assume one ancestral deme
        t_max = g.demes[0].epochs[0].end_time * 1.10
    t_min = min([d.epochs[-1].end_time for d in g.demes])
    if t_min == 0:
        t_min = 1
    if args.n is not None:
        n = args.n 
    else:
        n = 20
    if args.log_time:
        if t_min == 0:
            t_min = 1e4
        else:
            pass
        ts = np.logspace(np.log10(t_min), np.log10(t_max), n)
    else:
        ts = np.linspace(t_min, t_max, n)
    # augment with epoch boundaries
    epoch_edges = np.unique([e.start_time for d in g.demes for e in d.epochs][1:])
    ts = np.sort(np.concatenate((ts, epoch_edges)))
    rs = np.array([1e-6, 1e-5, 5e-5, 1e-3])
    u = args.u
    h2ts = {t: get_h2(g, t, u, rs, cross_pop=args.cross_pop) for t in ts}
    # construct H2 curves
    h2s = {}
    for t in h2ts: 
        for d in h2ts[t]:
            if d not in h2s:
                h2s[d] = np.full((len(rs), len(ts)), np.nan)
    for i, t in enumerate(h2ts):
        for d in h2ts[t]:
            h2s[d][:, i] = h2ts[t][d]
    print("DONE")
    fig, ax = plt.subplots(
        figsize=(6.5, 5.5), layout="constrained", sharey=True
    )
    colors = list(Category10_10) + list(TolRainbow23)
    linestyles = ["solid", "dashed", "dashdot", "dotted"]
    markers = ["^", "o", "d", "s"]
    ax.set_xlabel("$t$, years ago")
    ax.set_ylabel("$D^+(t)$")
    if args.data_file is not None:
        data = load_data(g, rs, args.data_file)
    else:
        data = None
    print(h2s.keys())
    for j, pop in enumerate(h2s):
        for i, r in enumerate(rs):
            if isinstance(pop, tuple):
                pop0, pop1 = pop
                label = f"{pop0},{pop1}"
                alpha = 0.5
            else:
                label = pop
                alpha = 1
            ax.plot(
                ts, h2s[pop][i], color=colors[j], linestyle=linestyles[i],
                label=label if i == 0 else None, alpha=alpha
            )
        if data is not None:
            if pop in data:
                for i, r in enumerate(rs):
                    ax.scatter(
                        [t_min], data[pop][i], color=colors[j], 
                        marker=markers[i], fc="none"
                    )
            elif isinstance(pop, tuple):
                if (pop[1], pop[0]) in data:
                    for i, r in enumerate(rs):
                        ax.scatter(
                            [t_min], data[(pop[1], pop[0]) ][i], color=colors[j], 
                            marker=markers[i], fc="none"
                        )
    ax.legend(fontsize=8, framealpha=0, ncols=2)
    handles = [
        mpl.lines.Line2D([], [], color="black", linestyle=linestyles[i],
        label=f"$r=${r}") for i, r in enumerate(rs)
    ]
    if args.data_file is not None:
        handles += [
            mpl.lines.Line2D([], [], color="black", marker=markers[j], 
            label=f"$r=${r}", linewidth=0, markerfacecolor="none") 
            for j, r in enumerate(rs)
        ]
        if args.log_time:
            ax.set_xscale("log")
            ax.set_xlim(t_min * 0.8,)
        else:
            ax.set_xlim(-0.05 * t_max,)
    else:
        if args.log_time:
            ax.set_xscale("log")
            ax.set_xlim(t_min,)
        else:
            ax.set_xlim(t_min,)
    if args.log_y:
        ax.set_yscale("log")
    fig.legend(handles=handles, framealpha=0, fontsize=8, ncols=len(rs),
               loc="lower center", bbox_to_anchor=(0.5, -0.10))
    if args.title is not None:
        fig.suptitle(args.title)
    ax.grid(alpha=0.3, axis="y")
    plt.savefig(args.out_file, dpi=244)

    return 


mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = \
    "\\usepackage{amsmath}\\usepackage{amssymb}"
mpl.rcParams['font.family'] = "serif"
mpl.rcParams['font.serif'] = "Computer Modern"
mpl.rcParams['savefig.bbox'] = "tight"
main()
