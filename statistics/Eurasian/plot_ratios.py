
import numpy as np
import pickle
import random
from dpluspy import bootstrapping, plotting


in_files = [f'sums/sums_chr{i}.pkl' for i in range(1, 23)]

regions = {}
for f in in_files:
    with open(f, 'rb') as fin:
        chrom_regions = pickle.load(fin)
    for key in chrom_regions:
        regions[key] = chrom_regions[key]
        pop_ids = chrom_regions[key]['pop_ids']
        bins = chrom_regions[key]['bins']


def weighted_bootstrap(
    regions,
    num_reps=None, 
    sample_size=None, 
    get_reps=False
):
    """
    Returns D+/H^2
    """
    if num_reps is None:
        num_reps = len(regions)
    if sample_size is None:
        sample_size = len(regions)
    means = bootstrapping.weighted_means_across_regions(regions)
    means[:-1] /= means[-1] ** 2

    labels = list(regions.keys())
    sample_size = len(regions)
    bootstrap_means = []
    for i in range(num_reps):
        samples = np.random.choice(labels, sample_size, replace=True)
        sampled_regions = [regions[sample] for sample in samples]
        _means = bootstrapping.weighted_means_across_replicates(sampled_regions)
        _means[:-1] /= _means[-1] ** 2
        bootstrap_means.append(_means)
    varcovs = []
    for i in range(len(means)):
        bin_means = np.array([_means[i] for _means in bootstrap_means])
        varcov_matrix = np.cov(bin_means.T)
        # this occurs when the bootstrap involves only one statistic
        if varcov_matrix.shape == ():
            varcov_matrix = varcov_matrix.reshape((1, 1))
        varcovs.append(varcov_matrix)
     
    return means, varcovs


means, varcovs = weighted_bootstrap(regions)
bins = np.loadtxt(bins)

#to_pops = ['Yoruba1', 'UstIshim', 'Stuttgart', 'Vindija', 'Denisova']
#means = bootstrapping.subset_means(means, pop_ids, to_pops=to_pops)
#varcovs = bootstrapping.subset_varcovs(varcovs, pop_ids, to_pops=to_pops)

plotting.plot_d_plus_curves(
    means=means,
    varcovs=varcovs,
    pop_ids=pop_ids,
    bins=bins,
    cols=11,
    aspect=1.2,
    ylabel='$D^+/\pi^2$',
    hline=1,
    out='fig_6_bootstrap_ratios.png',
)


