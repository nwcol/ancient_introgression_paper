
import pickle
from dpluspy import bootstrapping, plotting


stats = pickle.load(open('Eurasian_stats_weighted.pkl', 'rb'))
to_pops = stats['pop_ids'][:-2]
stats = bootstrapping.subset_statistics(stats, to_pops=to_pops, return_dict=True)

means = stats['means']
ratios = [m / means[-1] ** 2 for m in means[:-1]]
ratios.append(means[-1])
# This is not sound, I just want to have some measure of uncertainty on the plot
varcovs = stats['varcovs']
var_sqr_H = varcovs[-1] + means[-1] ** 2
ratio_varcovs = [c / var_sqr_H for c in varcovs[:-1]]
ratio_varcovs.append(varcovs[-1])

plotting.plot_d_plus_curves_comparison(
    means=ratios,
    varcovs=ratio_varcovs,
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    ax_size=2,
    cols=9,
    plot_grid=True,
    out='fig_5_ratios.png',
)


plotting.plot_d_plus_curves_comparison(
    means=ratios,
    varcovs=ratio_varcovs,
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    ax_size=2,
    cols=9,
    sharey=True,
    ylim=0,
    plot_grid=True,
    out='fig_3_ratios.png',
)

plotting.plot_d_plus_curves_comparison(
    means=stats['means'],
    varcovs=stats['varcovs'],
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    plot_H=True,
    ax_size=2,
    cols=9,
    out='fig_1_summary.png',
)

plotting.plot_d_plus_curves_comparison(
    means=stats['means'],
    varcovs=stats['varcovs'],
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    sharey=True,
    ax_size=2,
    cols=9,
    out='fig_2_summary.png',
)

