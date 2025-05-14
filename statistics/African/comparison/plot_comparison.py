
import pickle
from dpluspy import bootstrapping, plotting


pops = ['Yoruba1', 'KhomaniSan1', 'Mbuti1', 'Dinka1', 'Luhya1']

# These older statistics were parsed with different machinery and somewhat
# different masking, mutation weighting etc
_pops = ['Yoruba', 'Khomani_San', 'Mbuti', 'Dinka', 'Luhya']
older_stats = pickle.load(open('Africans_Behrer_10kb.pkl', 'rb'))['bootstrap']
_bins = older_stats['bins']
pop_ids = older_stats['pops']
means = [m for m in older_stats['means']]
varcovs = [m for m in older_stats['covs']]
means = bootstrapping.subset_means(means, pop_ids, to_pops=_pops)
varcovs = bootstrapping.subset_varcovs(varcovs, pop_ids, to_pops=_pops)

stats = pickle.load(open('../African_stats_weighted.pkl', 'rb'))
stats = bootstrapping.subset_statistics(stats, to_pops=pops, return_dict=True)
# stats_to_plot = utils._get_Dplus_names(pops)
# Hs_to_plot = utils._get_H_names(pops)
plotting.plot_d_plus_curves_comparison(
    means=[stats['means'], means],
    varcovs=[stats['varcovs'], varcovs],
    pop_ids=stats['pop_ids'],
    bins=[stats['bins'], _bins],
    plot_H=True,
    # stats_to_plot=stats_to_plot,
    # Hs_to_plot=Hs_to_plot,
    out='fig_1_summary.png',
)

plotting.plot_d_plus_curves_comparison(
    means=[stats['means'], means],
    varcovs=[stats['varcovs'], varcovs],
    pop_ids=stats['pop_ids'],
    bins=[stats['bins'], _bins],
    sharey=True,
    out='fig_2_summary.png',
)

