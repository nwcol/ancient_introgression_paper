
from bisect import bisect
import numpy as np
from h2py import util, h2_parsing


## two ideas. (1) revive site-pair u-weighting (2) bin site pairs by u prod


def bin_pair_products(
    mut_map, 
    rec_map, 
    bins, 
    left_lim=None,
    verbose=True
):
    if not left_lim:
        left_lim = len(rec_map)

    cum_vals = np.cumsum(mut_map)
    left_vals = mut_map[:left_lim]
    num_bins = len(bins) - 1
    products = np.zeros(num_bins, dtype=float)

    if bins[0] == 0:
        indices = np.arange(0, left_lim)
    else:
        indices = np.searchsorted(rec_map, rec_map[:left_lim] + bins[0]) - 1

    cum_sum0 = cum_vals[indices]

    for i, b in zip(range(num_bins), bins[1:]):
        indices = np.searchsorted(rec_map, rec_map[:left_lim] + b) - 1
        cum_sum1 = cum_vals[indices]
        cum_products = left_vals * (cum_sum1 - cum_sum0)
        products[i] = cum_products.sum()
        cum_sum0 = cum_sum1

        if verbose:
            print(util.get_time(), f"site pair products computed in bin {i}")

    return products


def double_binned(mut_map, rec_map, mut_bins, rec_bins, verbose=1000000):


    counts = np.zeros((len(rec_bins) - 1, len(mut_bins) - 1), dtype=np.float64)

    for i, (r_l, u_l) in enumerate(zip(rec_map, mut_map)):
        products = u_l * mut_map[i + 1:]
        distances = rec_map[i + 1:] - r_l
        counts += np.histogram2d(distances, products, bins=(rec_bins, mut_bins))

        if i % verbose == 0 and i > 0:
            print(util.get_time(), f"site pair products computed at locus {i}")

    return counts


def test_product_binning(mut_map, rec_map, bins, left_lim=None):
    ## naive and simple fxn for binning site mutation rate products by 
    ## recombination distance
    if left_lim is None:
        left_lim = len(rec_map)

    num_bins = len(bins) - 1
    binned_products = np.zeros(num_bins, dtype=np.float64)

    for i, (r_l, u_l) in enumerate(zip(rec_map[:left_lim], mut_map[:left_lim])):
        # we do not wish to count site i pairing with itself
        for r_r, u_r in zip(rec_map[i + 1:], mut_map[i + 1:]):
            distance = r_r - r_l
            bin_idx = bisect(bins, distance) - 1
            if bin_idx >= 0 and bin_idx < num_bins:
                binned_products[bin_idx] += u_l * u_r

    return binned_products
