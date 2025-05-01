## Functions for performing bootstraps and computing means.

from collections import defaultdict
import copy
import numpy as np
import pickle
import warnings

from . import utils


def load_raw_stats(filenames, load_mut_facs=False):
    """
    Load statistics from .pkl files. The expected format of each file is
    {
        "key1": {
            "sums": array([[...]]),
            "denoms": array([...])},
            "bins": array([...]),
            "pop_ids": [...], 
            "mut_facs": array([]) (optional)
            ...
        },
        "key2": {...},
        ...
    }

    :param filenames: Files from which to load raw statistics.

    :returns: population IDs, bin edges, and dictionary of raw data
    :rtypes: (list, np.ndarray, dict)
    """
    raw_data = defaultdict(dict)
    for filename in filenames:
        contents = pickle.load(open(filename, "rb+"))
        for key in contents:
            if key in raw_data: 
                raise ValueError(f"{key} appears twice in input")
            bins = contents[key]["bins"]
            pop_ids = contents[key]["pops"]
            raw_data[key]["sums"] = contents[key]["sums"]
            raw_data[key]["denoms"] = contents[key]["denoms"]
            raw_data[key]["weights"] = contents[key]["weights"]
            if load_mut_facs:
                if not "mut_facs" in contents[key]:
                    raise ValueError(f"Region {key} has no `mut_facs`")
                raw_data[key["mut_facs"]] = contents[key]["mut_facs"]

    return pop_ids, bins, raw_data


def load_bootstrap_means(filename, to_pops=None, size=None):
    # load list of bootstrap means
    with open(filename, "rb") as fin:
        contents = pickle.load(fin)
    if size is not None:
        all_means = contents["bootstrap_means"]
        samples = np.random.choice(
            np.arange(len(all_means)), size=size, replace=False
        )
        bootstrap_means = [all_means[i] for i in samples]
    else:
        bootstrap_means = contents["bootstrap_means"]
    pop_ids = contents["pop_ids"]
    if to_pops is None:
        ret_means = bootstrap_means
        ret_pop_ids = pop_ids
    else:
        ret_means = []
        for means in bootstrap_means:
            ret_means.append(utils.subset_means(means, pop_ids, to_pops))
        ret_pop_ids = to_pops
    bins = contents["bins"]

    return ret_means, bins, ret_pop_ids


def bootstrap(regions, num_reps=None, get_reps=False):
    """
    Perform a bootstrap to obtain covariance matrices for D+ and H statistics,
    estimated in genomic blocks. Operates upon sums of D+, H, and their
    respective denominators, so that regions are weighted appropriately. 
    
    :param regions: A dictionary with sums of D+ and H statistics from genomic
        regions as values, with arbitrary keys specifying region names.
    :type regions: dict
    :param num_reps: Number of bootstrap replicates to perform; if None, then 
        uses the number of regions (default None).
    :type num_reps: int
    :param get_reps: If True, return a list of bootstrap replicate means in 
        addition to means and covariances (default False).
    :type get_reps: bool
    
    :returns: Lists of mean and covariance arrays for each D+ bin and for the H
        statistics. Optionally also a list of bootstrap replicate means.
    :rtype: tuple (2 or 3-tuple of lists)
    """
    if num_reps is None:
        num_reps = len(regions)
    means = means_across_regions(regions)
    labels = list(regions.keys())
    sample_size = len(regions)
    bootstrap_means = []
    for i in range(num_reps):
        samples = np.random.choice(labels, sample_size, replace=True)
        sampled_regions = [regions[sample] for sample in samples]
        _means = means_across_replicates(sampled_regions)
        bootstrap_means.append(_means)
    varcovs = []
    for i in range(len(means)):
        bin_means = np.array([_means[i] for _means in bootstrap_means])
        varcov_matrix = np.cov(bin_means.T)
        # this occurs when the bootstrap involves only one statistic
        if varcov_matrix.shape == ():
            varcov_matrix = varcov_matrix.reshape((1, 1))
        varcovs.append(varcov_matrix)
    if get_reps:
        ret = (means, varcovs, bootstrap_means)
    else:
        ret = (means, varcovs)

    return ret


def means_across_regions(regions):
    """
    Compute mean statistics across genomic windows.
    
    :param regions: Dictionary of dictionaries that represent genomic windows,
        each containing summed D+, H statistics and respective denominators.
        Statistics should be numpy arrays, with the 0th dimension indexing bins.
    :type regions: dict

    :returns: A list holding the mean statistics in each bin.
    :rtype: list
    """
    sums = 0.0
    denoms = 0.0
    for key in regions:
        sums += regions[key]["sums"]
        denoms += regions[key]["denoms"]
    ext_denoms = np.repeat(denoms[:, np.newaxis], sums.shape[1], axis=1)
    raw_means = np.full(sums.shape, np.nan, dtype=np.float64)
    np.divide(sums, ext_denoms, where=ext_denoms > 0, out=raw_means)
    if np.any(np.isnan(raw_means)):
        warnings.warn("nan means exist in output")
    means = [raw_means[i] for i in range(len(raw_means))]

    return means


def means_across_replicates(replicates):
    """
    Compute mean statistics across a list of replicates. 
    
    :param replicates: List of dictionaries that hold summed D+, H statistics
        and the respective denominators as numpy arrays.
    :type replicates: list

    :returns: A list of mean statistics for each bin.
    :rtype: list
    """
    rep_dict = {i: replicate for i, replicate in enumerate(replicates)}
    means = means_across_regions(rep_dict)

    return means



def weighted_bootstrap(
    regions,
    num_reps=None, 
    sample_size=None, 
    get_reps=False
):
    """
    Perform a bootstrap on a dictionary of region-specific D+ and H sums, 
    using the mutation rate-weighted estimator to de-distort the shape of the
    D+ curve. 

    Each region should be represented by a weighted dictionary in `regions`,
    minimally containing 'sums', 'num_sites' and 'mut_facs'. 
    """
    if num_reps is None:
        num_reps = len(regions)
    if sample_size is None:
        sample_size = len(regions)
    means = weighted_means_across_regions(regions)
    labels = list(regions.keys())
    sample_size = len(regions)
    bootstrap_means = []
    for i in range(num_reps):
        samples = np.random.choice(labels, sample_size, replace=True)
        sampled_regions = [regions[sample] for sample in samples]
        _means = weighted_means_across_replicates(sampled_regions)
        bootstrap_means.append(_means)
    varcovs = []
    for i in range(len(means)):
        bin_means = np.array([_means[i] for _means in bootstrap_means])
        varcov_matrix = np.cov(bin_means.T)
        # this occurs when the bootstrap involves only one statistic
        if varcov_matrix.shape == ():
            varcov_matrix = varcov_matrix.reshape((1, 1))
        varcovs.append(varcov_matrix)
    if get_reps:
        ret = (means, varcovs, bootstrap_means)
    else:
        ret = (means, varcovs)

    return ret


def weighted_means_across_regions(regions):
    """
    Compute mutation-rate weighted D+ across a dictionary of regions.
    """
    sums = 0.0
    pair_counts = 0.0
    num_sites = 0.0
    mut_prods = 0.0
    for key in regions:
        sums += regions[key]['sums']
        pair_counts += regions[key]['denoms'][:-1]
        num_sites += regions[key]['denoms'][-1]
        mut_prods += regions[key]['mut_facs'][:-1]
    # Compute the u-weighted denominator
    facs = mut_prods / (mut_prods.sum() / pair_counts.sum())
    weighted_denoms = mut_prods / facs
    denoms = np.append(weighted_denoms, num_sites)
    
    ext_denoms = np.repeat(denoms[:, np.newaxis], sums.shape[1], axis=1)
    raw_means = np.full(sums.shape, np.nan, dtype=np.float64)
    np.divide(sums, ext_denoms, where=ext_denoms > 0, out=raw_means)
    if np.any(np.isnan(raw_means)):
        warnings.warn("nan means exist in output")
    means = [raw_means[i] for i in range(len(raw_means))]

    return means


def weighted_means_across_replicates(replicates):
    """
    Operates on a list of dictionaries; wraps `weighted_means_across_regions`.
    """
    rep_dict = {i: replicate for i, replicate in enumerate(replicates)}
    means = weighted_means_across_regions(rep_dict)

    return means


### DEPRECATED. this stuff doesn't work well!!


def compute_weighted_denoms(mut_facs, denoms):
    """
    Compute weighted denominators for D+.

    Weighting takes the form of adjusting the effective number of locus pairs
    (the denominator) of each bin in inverse proportion to its mean product
    of locus mutation rates u_L * u_R (relative to the average across bins).

    :param mut_facs: Binned sums of locus-pair mutation rates u_L * u_R. The 
        last element should be the sum of locus mutation rates, which is not 
        used here.
    :param denoms: Binned counts of locus pairs. It is assumed that the last
        element contains the denominator for H (the number of loci), which is
        ignored here and returned as the last element of `weights`.

    :returns: Weights for adjusting D+
    :rtype: np.ndarray
    """
    num_sites = denoms[-1]
    mut_prods = mut_facs[:-1]
    locus_pairs = denoms[:-1]
    #factor = mut_prods.sum() / locus_pairs.sum()
    factor = (mut_facs[-1] / num_sites) ** 2
    _weighted_denoms = mut_prods / factor 
    weighted_denoms = np.append(_weighted_denoms, num_sites)

    return weighted_denoms


def compute_weights(mut_facs, denoms):
    """
    Compute the ratio of the weighted D+ denominator over the unweighted one.
    """
    weighted_denoms = compute_weighted_denoms(mut_facs, denoms)
    weights = weighted_denoms / denoms

    return weights



