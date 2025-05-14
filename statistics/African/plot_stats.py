
import pickle
from dpluspy import bootstrapping, plotting


stats = pickle.load(open('African_stats_weighted.pkl', 'rb'))
pops = ['Yoruba1', 'Yoruba2', 'KhomaniSan1', 'Mbuti1', 'Dinka1', 'Luhya1']
stats = bootstrapping.subset_statistics(stats, to_pops=pops, return_dict=True)
# stats_to_plot = utils._get_Dplus_names(pops)
# Hs_to_plot = utils._get_H_names(pops)

means = stats['means']
ratios = [m / means[-1] ** 2 for m in means[:-1]]
ratios.append(means[-1])
# This is not sound, I just want to have some measure of uncertainty on the plot
varcovs = stats['varcovs']
var_sqr_H = varcovs[-1] + means[-1] ** 2
ratio_varcovs = [c / var_sqr_H for c in varcovs[:-1]]
ratio_varcovs.append(varcovs[-1])

plotting.plot_d_plus_curves_comparison(
    means=stats['means'],
    varcovs=stats['varcovs'],
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    sharey=True,
    plot_grid=True,
    cols=6,
    out='fig_2_summary.png',
)
          
plotting.plot_d_plus_curves_comparison(
    means=ratios,
    varcovs=ratio_varcovs,
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    cols=6,
    sharey=True,
    ylim=0,
    plot_grid=True,
    out='fig_4_ratios.png',
)

plotting.plot_d_plus_curves_comparison(
    means=stats['means'],
    varcovs=stats['varcovs'],
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    plot_grid=True,
    cols=6,
    out='fig_3_summary.png',
)

plotting.plot_d_plus_curves_comparison(
    means=stats['means'],
    varcovs=stats['varcovs'],
    pop_ids=stats['pop_ids'],
    bins=stats['bins'],
    plot_H=True,
    cols=6,
    # stats_to_plot=stats_to_plot,
    # Hs_to_plot=Hs_to_plot,
    out='fig_1_summary.png',
)


