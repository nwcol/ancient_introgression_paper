
import numpy as np
import dpluspy


def compute_numer_onepop(
    genotypes, 
    rec_map,
    bins, 
    mut_map, 
    u_bar, 
    verbose=False
):
    """
    Compute the numerator of the u-adjusted statistic. This is

    u_bar ** 2 * sum_(i,j) D+(i,j) / (u_i * u_j),

    and the whole estimator is

    u_bar ** 2 / n_pairs * sum_(i,j) D+(i,j) / (u_i * u_j)

    :param genotypes: Array of genotypes for a single diploid
    :param rec_map: Array of recombination map coordinates for genotyped sites
    :param bins: Array of recombination bin edges
    :param mut_map: Array of estimated mutation rates at genotyped sites
    """
    # indicator for one-locus heterozygosity
    indicator = 1.0 * (genotypes[:, 0] != genotypes[:, 1])
    weights = u_bar / mut_map
    indicator *= weights

    stats = dpluspy.parsing._count_locus_pairs(
        rec_map, 
        bins,
        weights=indicator,
        verbose=verbose
    )

    return stats


def compute_numer_crosspop():


    return
