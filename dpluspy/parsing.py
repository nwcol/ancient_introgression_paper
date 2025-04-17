## Houses functions for estimating D+ from sequence data.

from collections import defaultdict
import copy
import gzip
import numpy as np
import re
import scipy
import warnings

from . import utils


def compute_statistics(
    vcf_file,
    bed_file=None,
    pop_file=None,
    pop_mapping=None,
    rec_map_file=None,
    L=None,
    r=None,
    map_type="linear",
    interval=None,
    interval_between=None,
    r_bins=None,
    phased=False,
    cross_pop=True
):
    """
    Compute the D+ and H statistics from a .vcf file. 

    :param interval_between: 2-tuple or list of 2-tuples or lists, specifying
        the interval of the left and right loci when parsing D+ between two
        windows. H is left as 0 in this case (default None).
    """
    if interval is not None:
        within = True
        assert len(interval) == 2
    elif interval_between is not None:
        within = False
        assert len(interval_between) == 2
        left_interval = interval_between[0]
        right_interval = interval_between[1]
        assert len(left_interval) == len(right_interval) == 2
        assert right_interval[0] >= left_interval[1]
        interval = (left_interval[0], right_interval[1])
    else:
        within = True
        warnings.warn("no interval given; parsing all sites")

    if pop_file is not None and pop_mapping is not None:
        raise ValueError('You cannot use both `pop_file` and `pop_mapping`')
    if pop_file is not None:
        pop_mapping = load_pop_file(pop_file)

    if r is not None:
        if L is None:
            if interval is not None:
                L = interval[-1]
            else:
                raise ValueError("please provide L")
        rec_map = get_uniform_recombination_map(L, r, kind=map_type)
    else:
        rec_map = load_recombination_map(rec_map_file, kind=map_type)

    if isinstance(r_bins, str):
        r_bins = np.loadtxt(r_bins)
    bins = utils.map_function(r_bins)

    if within:
        sites, genotype_matrix, sample_ids = read_genotypes(
            vcf_file, bed_file=bed_file, interval=interval
        )
        pop_genotypes = build_pop_genotypes(
            genotype_matrix, sample_ids, pop_mapping=pop_mapping
        )
        site_map = rec_map(sites)
        sums = compute_stats_within(
            pop_genotypes, site_map, bins, cross_pop=cross_pop, phased=phased,
        )
        pop_ids = list(pop_genotypes.keys())
    else:
        sites_left, genotype_matrix_left, sample_ids = read_genotypes(
            vcf_file, bed_file=bed_file, interval=left_interval
        )
        pop_genotypes_left = build_pop_genotypes(
            genotype_matrix_left, sample_ids, pop_mapping=pop_mapping
        )
        sites_right, genotype_matrix_right, _ = read_genotypes(
            vcf_file, bed_file=bed_file, interval=right_interval
        )
        pop_genotypes_right = build_pop_genotypes(
            genotype_matrix_right, sample_ids, pop_mapping=pop_mapping
        )
        site_map_left = rec_map(sites_left)
        site_map_right = rec_map(sites_right)
        sums = compute_stats_between(
            pop_genotypes_left, pop_genotypes_right, 
            site_map_left, site_map_right,
            bins, cross_pop=True, phased=False
        )
        pop_ids = list(pop_genotypes_left.keys())

    return sums


def compute_denominators(
    bed_file=None,
    rec_map_file=None,
    map_type="linear",
    mut_map_file=None,
    mut_map_col='u',
    L=None,
    r=None,
    interval=None,
    interval_between=None,
    r_bins=None
):
    """
    Compute denominators for the D+ and H statistics. 


    """
    if interval is not None:
        within = True
        assert len(interval) == 2
    elif interval_between is not None:
        within = False
        assert len(interval_between) == 2
        left_interval = interval_between[0]
        right_interval = interval_between[1]
        assert len(left_interval) == len(right_interval) == 2
        assert right_interval[0] >= left_interval[1]
        interval = (left_interval[0], right_interval[1])
    else:
        within = True
        warnings.warn("no interval given; parsing all sites")

    if r is not None:
        if L is None:
            if interval is not None:
                L = interval[-1]
            else:
                raise ValueError("please provide L")
        rec_map = get_uniform_recombination_map(L, r, kind=map_type)
    else:
        rec_map = load_recombination_map(rec_map_file, kind=map_type)

    if bed_file is not None:
        all_positions = utils.read_bedfile_positions(bed_file)
    else:
        if L is None:
            raise ValueError("L or a bed_file is required to compute denoms")
        all_positions = np.arange(L)

    if isinstance(r_bins, str):
        r_bins = np.loadtxt(r_bins)
    bins = utils.map_function(r_bins)

    if within:
        if interval is None:
            positions = all_positions
        else:
            positions = all_positions[
                (all_positions >= interval[0]) & (all_positions < interval[1])
            ]
        pos_map = rec_map(positions)
        denoms = count_locus_pairs(pos_map, bins, verbose=True)
        denoms = np.append(denoms, len(positions))

        if mut_map_file is not None:
            if mut_map_file.endswith(".npy"):
                _mut_map = np.load(mut_map_file)
                mut_map = _mut_map[positions]
            elif ".bedgraph" in mut_map_file:
                windows, mut_data = utils.read_bedgraph(mut_map_file)
                avg_mut = mut_data[mut_map_col]
                mut_map = mut_map_discretized(windows, avg_mut, positions)
            mut_facs= count_locus_pairs(
                pos_map, bins, weights=mut_map, verbose=True
            )
            sum_mut = np.sum(mut_map)
            mut_facs = np.append(mut_facs, sum_mut)
        else:
            mut_facs = None
    else:
        left_positions = all_positions[
            (all_positions >= left_interval[0]) 
            & (all_positions < left_interval[1])
        ]
        pos_map_left = rec_map(left_positions)
        right_positions = all_positions[
            (all_positions >= right_interval[0]) 
            & (all_positions < right_interval[1])
        ]
        pos_map_right = rec_map(right_positions)
        denoms = count_locus_pairs_between(
            pos_map_left, pos_map_right, bins, verbose=True
        )
        denoms = np.append(denoms, 0)
        if mut_map_file is not None:
            if mut_map_file.endswith(".npy"):
                _mut_map = np.load(mut_map_file)
                mut_map_left = _mut_map[left_positions]
                mut_map_right = _mut_map[right_positions]
            elif ".bedgraph" in mut_map_file:
                windows, mut_data = utils.read_bedgraph(mut_map_file)
                avg_mut = mut_data[mut_map_col]
                mut_map_left = mut_map_discretized(
                    windows, avg_mut, left_positions
                )
                mut_map_right = mut_map_discretized(
                    windows, avg_mut, right_positions
                )
            mut_facs = count_locus_pairs_between(
                pos_map_left, pos_map_right, bins, 
                weights_l=mut_map_left, weights_r=mut_map_right
            )
            mut_facs = np.append(mut_facs, 0)
        else:
            mut_facs = None 

    if mut_facs is None:
        return denoms
    else:
        return denoms, mut_facs
    

def load_pop_file(pop_file):
    """
    Load a population file.
    """
    pop_mapping = defaultdict(list)
    with open(pop_file, 'r') as fin:
        for line in fin:
            sample, pop = line.split()
            pop_mapping[pop].append(sample)
    return pop_mapping


def build_pop_genotypes(gt_matrix, sample_ids, pop_mapping=None):
    """
    From an array of genotypes encoding 
    
    :param pop_mapping: Dictionary, mapping population IDs to lists of sample
        IDs. 
    """
    if pop_mapping is None:
        pop_mapping = {sample_id: [sample_id] for sample_id in sample_ids}
    pop_ids = list(pop_mapping.keys())
    pop_indices = {}
    for pop_id in pop_ids:
        samples = pop_mapping[pop_id]
        pop_indices[pop_id] = [sample_ids.index(sample) for sample in samples]
    pop_genotypes = {}
    for pop_id in pop_ids:
        pop_genotypes[pop_id] = gt_matrix[:, pop_indices[pop_id]]

    return pop_genotypes 


def flatten_genotypes(pop_genotypes):
    """
    Convert a dictionary of arrays with shapes (s, n, 2) to a dictionary of
    arrays with shape (s, 2 * n). These may have the interpretation of being
    an array of haplotypes, but in some contexts it is also just convenient to
    have genotypes represented in this way (e.g. estimating pairwise diversity)
    """
    flat_genotypes = {}
    for pop_id in pop_genotypes:
        array = pop_genotypes[pop_id]
        s, n, _ = array.shape
        flat_array = np.reshape(array, (s, 2 * n))
        flat_genotypes[pop_id] = flat_array

    return flat_genotypes


def compute_pi(pop_genotypes, cross_pop=True):
    """
    Compute nucleotide diversity in a contiguous genomic block. Returns an 
    array of sums (to be normalized by L).
    """
    flat_genotypes = flatten_genotypes(pop_genotypes)
    pop_ids = list(flat_genotypes.keys())
    num_pops = len(pop_ids)
    if cross_pop:
        num_stats = (num_pops + num_pops ** 2) // 2
    else:
        num_stats = num_pops
    sums = np.zeros(num_stats, dtype=np.float64)
    idx = 0
    for i, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[i:]:
            if pop_i == pop_j: 
                alleles = flat_genotypes[pop_i]
                _, n = alleles.shape
                numer = 0.0
                for k in range(n - 1):
                    for l in range(k + 1, n):
                        numer += (alleles[:, k] != alleles[:, l]).sum()
                sum_i = numer / (n * (n - 1) / 2)
                sums[idx] = sum_i
            else:
                if not cross_pop:
                    continue
                alleles_i = flat_genotypes[pop_i]
                alleles_j = flat_genotypes[pop_j]
                _, ni = alleles_i.shape
                _, nj = alleles_j.shape
                numer = 0.0
                for k in range(ni):
                    for l in range(nj):
                        numer += (alleles_i[:, k] != alleles_j[:, l]).sum()
                sum_ij = numer / (ni * nj)
                sums[idx] = sum_ij 
            idx += 1

    return sums


def compute_stats_within(
    pop_genotypes, 
    site_map, 
    bins, 
    cross_pop=True,
    phased=False,
):
    """
    Compute statistics in a contiguous genomic block.

    :param pop_genotypes: Should instead hold haplotypes if `phased` is True.
    """
    pop_ids = list(pop_genotypes.keys())
    num_pops = len(pop_ids)
    if cross_pop:
        num_stats = (num_pops + num_pops ** 2) // 2
    else:
        num_stats = num_pops
    sums = np.zeros((len(bins), num_stats))
    idx = 0
    for i, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[i:]:
            if pop_i == pop_j:
                Gt = pop_genotypes[pop_i]
                if phased:
                    sums[:-1, idx] = haplotype_Dplus(Gt, site_map, bins)
                else:
                    sums[:-1, idx] = genotype_Dplus(Gt, site_map, bins)
            else:
                if not cross_pop:
                    continue
                Gi = pop_genotypes[pop_i]
                Gj = pop_genotypes[pop_j]
                if phased:
                    sums[:-1, idx] = cross_haplotype_Dplus(Gi, Gj, site_map, bins)
                else:
                    sums[:-1, idx] = cross_genotype_Dplus(Gi, Gj, site_map, bins)
            idx += 1
    sums[-1] = compute_pi(pop_genotypes, cross_pop=cross_pop)

    return sums


def compute_stats_between(
    pop_genotypes_l,
    pop_genotypes_r,
    site_map_l,
    site_map_r,
    bins, 
    cross_pop=True,
    phased=False
):
    """
    Compute statistics between two contiguous genomic blocks.
    """
    pop_ids = list(pop_genotypes_l.keys())
    num_pops = len(pop_ids)
    if cross_pop:
        num_stats = (num_pops + num_pops ** 2) // 2
    else:
        num_stats = num_pops
    sums = np.zeros((len(bins), num_stats))
    idx = 0
    for i, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[i:]:
            if pop_i == pop_j:
                G_l = pop_genotypes_l[pop_i]
                G_r = pop_genotypes_r[pop_i]
                if phased:
                    sums[:-1, idx] = haplotype_Dplus_between(
                        G_l, G_r, site_map_l, site_map_r, bins
                    )
                else:
                    sums[:-1, idx] = genotype_Dplus_between(
                        G_l, G_r, site_map_l, site_map_r, bins
                    )
            else:
                if not cross_pop:
                    continue
                G_li = pop_genotypes_l[pop_i]
                G_lj = pop_genotypes_l[pop_j]
                G_ri = pop_genotypes_r[pop_i]
                G_rj = pop_genotypes_r[pop_j]
                if phased:
                    sums[:-1, idx] = cross_haplotype_Dplus_between(
                        G_li, G_lj, G_ri, G_rj, site_map_l, site_map_r, bins
                    )
                else:
                    sums[:-1, idx] = cross_genotype_Dplus_between(
                        G_li, G_lj, G_ri, G_rj, site_map_l, site_map_r, bins
                    )
            idx += 1
    sums[-1] = 0

    return sums


## Estimators


def haplotype_Dplus(haplotypes, site_map, bins):
    """
    
    """
    n = haplotypes.shape[1]
    if n == 2:
        weights = haplotypes[:, 0] != haplotypes[:, 1]
        Dplus = count_locus_pairs(site_map, bins, weights=weights)
    else:
        numer = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                numer += haplotype_Dplus(haplotypes[:, [i, j]], site_map, bins)
        Dplus = numer / (n * (n - 1) / 2)

    return Dplus


def haplotype_Dplus_between(
    haplotypes_l, 
    haplotypes_r, 
    site_map_l, 
    site_map_r, 
    bins
):
    """
    
    """
    n = haplotypes_l.shape[1]
    if n == 2:
        weights_l = haplotypes_l[:, 0] != haplotypes_l[:, 1]
        weights_r = haplotypes_r[:, 0] != haplotypes_r[:, 1]
        Dplus = count_locus_pairs_between(
            site_map_l,
            site_map_r, 
            bins, 
            weights_l=weights_l,
            weights_r=weights_r
        )
    else:
        numer = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                numer += haplotype_Dplus_between(
                    haplotypes_l[:, [i]], 
                    haplotypes_r[:, [j]], 
                    site_map_l,
                    site_map_r, 
                    bins
                )
        Dplus = numer / (n * (n - 1) / 2)

    return Dplus


def cross_haplotype_Dplus(haplotypes_i, haplotypes_j, site_map, bins):
    """
    
    """
    ni = haplotypes_i.shape[1]
    nj = haplotypes_j.shape[1]
    if ni == 1 and nj == 1:
        weights = haplotypes_i[:, 0] != haplotypes_j[:, 0]
        Dplus = count_locus_pairs(site_map, bins, weights=weights)
    else:
        numer = 0.0
        for k in range(ni):
            for l in range(nj):
                numer += cross_haplotype_Dplus(
                    haplotypes_i[:, [k]], haplotypes_j[:, [l]], site_map, bins
                )
        Dplus = numer / (ni * nj)

    return Dplus
 

def cross_haplotype_Dplus_between(
    haplotypes_li, 
    haplotypes_lj,
    haplotypes_ri, 
    haplotypes_rj,
    site_map_l, 
    site_map_r, 
    bins
):
    """
    
    """
    ni = haplotypes_li.shape[1]
    nj = haplotypes_lj.shape[1]
    if ni == 1 and nj == 1:
        weights_l = haplotypes_li[:, 0] != haplotypes_lj[:, 1]
        weights_r = haplotypes_ri[:, 0] != haplotypes_rj[:, 1]
        Dplus = count_locus_pairs_between(
            site_map_l,
            site_map_r, 
            bins, 
            weights_l=weights_l,
            weights_r=weights_r
        )
    else:
        numer = 0.0
        for k in range(ni):
            for l in range(nj):
                numer += haplotype_Dplus_between(
                    haplotypes_li[:, [k]], 
                    haplotypes_lj[:, [l]],
                    haplotypes_ri[:, [k]], 
                    haplotypes_rj[:, [l]],
                    site_map_l,
                    site_map_r, 
                    bins
                )
        Dplus = numer / (ni * nj)

    return Dplus


def genotype_Dplus(genotypes, site_map, bins):
    """
    Supports sites with >2 alleles.

    :param genotypes: Array with shape (s, n, 2)  
    """
    n = genotypes.shape[1]
    if n == 1:
        weights = genotypes[:, 0, 0] != genotypes[:, 0, 1]
        Dplus = count_locus_pairs(site_map, bins, weights=weights)
    else:
        numer = 0.0
        for i in range(n):
            numer += genotype_Dplus(genotypes[:, [i], :], site_map, bins)
        Dplus = numer / n

    return Dplus


def genotype_Dplus_between(
    genotypes_l,
    genotypes_r,
    site_map_l,
    site_map_r,
    bins
):
    """
    Compute binned D+ sums between two genomic blocks. 


    """
    n = genotypes_l.shape[1]
    if n == 1:
        weights_l = genotypes_l[:, 0, 0] != genotypes_l[:, 0, 1]
        weights_r = genotypes_r[:, 0, 0] != genotypes_r[:, 0, 1]
        Dplus = count_locus_pairs_between(
            site_map_l, 
            site_map_r, 
            bins,
            weights_l=weights_l, 
            weights_r=weights_r
        )
    else:
        numer = 0.0
        for i in range(n):
            numer += genotype_Dplus_between(
                genotypes_l[:, [i], :], 
                genotypes_r[:, [i], :],
                site_map_l,
                site_map_r,
                bins
            )
        Dplus = numer / n

    return Dplus


def cross_genotype_Dplus(genotypes_i, genotypes_j, site_map, bins):
    """
    
    """
    ni = genotypes_i.shape[1]
    nj = genotypes_j.shape[1]
    if ni == 1 and nj == 1:
        weights = pi_xy(genotypes_i[:, 0], genotypes_j[:, 0])
        Dplus = count_locus_pairs(site_map, bins, weights=weights)
    else:
        numer = 0.0
        for k in range(ni):
            for l in range(nj):
                numer += cross_genotype_Dplus(
                    genotypes_i[:, [k], :], 
                    genotypes_j[:, [l], :], 
                    site_map, 
                    bins
                )
        Dplus = numer / (ni * nj)

    return Dplus


def cross_genotype_Dplus_between(
    genotypes_li,
    genotypes_lj,
    genotypes_ri,
    genotypes_rj,
    site_map_l,
    site_map_r,
    bins
):
    """
    
    """
    ni = genotypes_li.shape[1]
    nj = genotypes_lj.shape[1]
    if ni == 1 and nj == 1:
        weights_l = pi_xy(genotypes_li[:, 0], genotypes_lj[:, 0])
        weights_r = pi_xy(genotypes_ri[:, 0], genotypes_rj[:, 0])
        Dplus = count_locus_pairs_between(
            site_map_l, 
            site_map_r, 
            bins,
            weights_l=weights_l, 
            weights_r=weights_r
        )
    else:
        numer = 0.0
        for k in range(ni):
            for l in range(nj):
                numer += cross_genotype_Dplus(
                    genotypes_li[:, [k], :],
                    genotypes_lj[:, [l], :],
                    genotypes_ri[:, [k], :],
                    genotypes_rj[:, [l], :],
                    site_map_l,
                    site_map_r,
                    bins
                )
        Dplus = numer / (ni * nj)

    return Dplus


def pi_xy(genotypes_i, genotypes_j):
    """
    Compute the pairwise divergence between two diploids. This is the nucleotide 
    diversity, conditional on sampling one allele copy from each diploid.

    :param genotypes_i: Array of allelic states with shape (s, 2) for diploid i.
    :param genotypes_j: Array of allelic states for diploid j with shape (s, 2)

    :returns: Array of pi_ij with shape (s,)
    :rtype: np.ndarray 
    """
    pairwise_diff = genotypes_i[:, :, np.newaxis] != genotypes_j[:, np.newaxis]
    pi = pairwise_diff.sum((2, 1)) / 4

    return pi


## Computing denominators and weights


def compute_mut_facs(pos_map, bins, mut_map):

    mut_prods = count_locus_pairs(pos_map, bins, weights=mut_map, verbose=True)
    sum_mut = np.sum(mut_map)
    mut_facs = np.append(mut_prods, sum_mut)

    return mut_facs


def mut_map_discretized(intervals, rates, positions):
    """
    Assign mutation rates to `positions`
    """
    site_rates = np.zeros(intervals[-1, 1])
    for rate, (start, end) in zip(rates, intervals):
        site_rates[start:end] = rate
    subset_site_rates = site_rates[positions]

    return subset_site_rates


## Locus pair-counting functions


def count_locus_pairs(site_map, bins, weights=None, verbose=False):
    """
    Computes the numbers of site pairs that fall within each of a series of 
    recombination bins, in a contiguous genomic region. 

    Used to compute D+ and its denominator. 

    :param site_map: Array giving the recombination map coordinates of sites
        in linear units (cM or M).
    :param bins: Array of recombination bin edges, given in the same unit as 
        the map (cM or M). 
    :weights: An array with length equal to `sitemap` assigning a weight to each 
        site (default None). Computing counts without weights is equivalent to 
        giving every site weight 1.
    :returns: Array of binned locus pair counts.
    """
    num_bins = len(bins) - 1

    if weights is not None:
        if bins[0] == 0:
            indices = np.arange(1, len(site_map) + 1)
        else:
            indices = np.searchsorted(site_map, site_map + bins[0])
        cum_weights = np.concatenate(([0], np.cumsum(weights)))
        cum_sum0 = cum_weights[indices]
        sums = np.zeros(num_bins, dtype=np.float64)
        for i, b in enumerate(bins[1:]):
            indices = np.searchsorted(site_map, site_map + b)
            cum_sum1 = cum_weights[indices]
            sums[i] = (weights * (cum_sum1 - cum_sum0)).sum()
            cum_sum0 = cum_sum1
            if verbose:
                print(utils.get_time(), f"locus pairs summed (within) in bin {i}")

    # TODO add a thing that returns 0s when all pair distances exceed highest bin edge

    else:
        if bins[0] == 0:
            edge0 = np.arange(1, len(site_map) + 1)
        else:
            edge0 = np.searchsorted(site_map, site_map + bins[0])
            assert np.all(edge0 > 0)
        sums = np.zeros(num_bins, dtype=np.int64)
        for i, b in enumerate(bins[1:]):
            edge1 = np.searchsorted(site_map, site_map + b)
            sums[i] = (edge1 - edge0).sum() 
            edge0 = edge1
            if verbose:
                print(utils.get_time(), f"locus pairs summed (within) in bin {i}")

    return sums


def count_locus_pairs_between(
    site_map_l, 
    site_map_r, 
    bins, 
    weights_l=None, 
    weights_r=None,
    verbose=False
):
    """
    Computes the numbers of site pairs that fall within each of a series of 
    recombination bins between two separate and internally contiguous genomic
    blocks. The left block `1` must have lower map coordinates than the right
    block. 

    Used to compute D+ and its denominator. 

    :params sitemap1, sitemap2: Array giving the recombination map coordinates 
        of sites in linear units (cM or M) for the left and right blocks.
    :param bins: Array of recombination bin edges, given in the same unit as 
        the map (cM or M). 
    :params weights1, weights1: Array with lengths equal to `sitemap1` and 
        `sitemap2` respectively, assigning weight site in each block 
        (default None).
    :returns: Array of binned locus pair counts.
    """
    if not np.max(site_map_l) <= np.min(site_map_r):
        raise ValueError("Block 1 have lower coordinate than block 2")
    if (weights_l is not None) ^ (weights_r is not None):
        raise ValueError("You must provide weights for both blocks")
    if weights_l is not None:
        if len(weights_l) != len(site_map_l):
            raise ValueError("Map and weight lengths mismatch for block 1")
        if len(weights_r) != len(site_map_r):
            raise ValueError("Map and weight lengths mismatch for block 2")

    num_bins = len(bins) - 1

    # TODO add a thing that returns 0s when all pair distances exceed highest bin edge

    if weights_l is not None:
        indices = np.searchsorted(site_map_r, site_map_l + bins[0])
        assert np.all(indices >= 0)
        cum_weights2 = np.concatenate(([0], np.cumsum(weights_r)))
        cum_sum0 = cum_weights2[indices]
        sums = np.zeros(num_bins, dtype=np.float64)
        for i, b in enumerate(bins[1:]):
            indices = np.searchsorted(site_map_r, site_map_l + b)
            assert np.all(indices >= 0)
            cum_sum1 = cum_weights2[indices]
            sums[i] = (weights_l * (cum_sum1 - cum_sum0)).sum()
            cum_sum0 = cum_sum1
            if verbose:
                print(utils.get_time(), f"locus pairs summed (between) in bin {i}")

    else:
        edge0 = np.searchsorted(site_map_r, site_map_l + bins[0])
        sums = np.zeros(num_bins, dtype=np.int64)
        for i, b in enumerate(bins[1:]):
            edge1 = np.searchsorted(site_map_r, site_map_l + b)
            sums[i] = (edge1 - edge0).sum() 
            edge0 = edge1
            if verbose:
                print(utils.get_time(), f"locus pairs summed (between) in bin {i}")

    return sums


## Utilities


def get_uniform_recombination_map(L, r, kind="linear"):
    """
    Generate a function that interpolates map coordinates for a uniform 
    recombination with rate `r` and length `L`. 

    :param L: Length of the map.
    :param r: Map rate, in units of r (recombination frequency).
    :param kind: The type of interpolation to use (default 'linear').

    :returns: Function that interpolates for a uniform map
    :rtype: scipy.interpolate.interp1d 
    """
    if kind not in ("nearest", "linear", "previous", "next"):
        raise ValueError("Invalid `kind`")
    x = np.arange(L)
    y = r * 100 * x
    mapfunc = scipy.interpolate.interp1d(
        x, y, kind=kind, bounds_error=False, fill_value=(y[0], y[-1])
    )
    return mapfunc


def load_recombination_map(file, map_col="Map(cM)", kind="linear"):
    """
    Load a recombination map and return an interpolate function.

    :param file: Filename of recombination map.
    :param map_col: Title of column containing map coordinates.
    :param kind: The type of interpolation to use (default 'linear').

    :returns: Interpolate function
    :rtype: scipy.interpolate.interp1d 
    """
    if kind not in ("nearest", "linear", "previous", "next"):
        raise ValueError("Invalid `kind`")
    x, y = utils.read_hapmap_rec_map(file, map_col=map_col)
    mapfunc = scipy.interpolate.interp1d(
        x, y, kind=kind, bounds_error=False, fill_value=(y[0], y[-1])
    )
    return mapfunc


def read_genotypes(
    vcf_file, 
    bed_file=None, 
    multiallelic=False,
    missing_to_ref=True,
    interval=None
):
    """
    Read sites and genotypes from a .vcf file.

    If return_dict is True, returns a dictionary mapping sites to site genotype
    arrays with shapes (n, 2) where n is the number of samples. Otherwise,
    returns an array of sites, an array of genotypes with shape (s, n, 2) where
    s is the number of sites. 

    We encode genotypes represented as A1/A2 or A1|A2 in a .vcf in the form 
    [A1, A2]. Thus if data is phased, the genotype array can be converted into
    a haplotype array by flattening the dimension.

    :param vcf_file: Filename of a .vcf file.
    :param bed_file: Filename for .bed filter to impose on sites (default None).
    :param multiallelic: If True, do not skip multiallelic sites (default False)
    :param missing_to_ref: If True, genotypes ./. and .|. will be read as 0/0
        or 0|0 respectively. Default False skips sites with any missing data.
    :param interval: 2-tuple or list specifying upper and lower bounds on 
        positions (default None).

    :returns: Array of 0-indexed sites, array of genotypes, list of sample IDs
    """
    if bed_file is not None:
        regions = utils.read_bedfile(bed_file)
        mask = utils.regions_to_mask(regions)
        masked = True
    else:
        masked = False

    open_func = gzip.open if vcf_file.endswith('.gz') else open
    _sites = []
    _genotypes = []

    with open_func(vcf_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode().strip()
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    sample_ids = line.split()[9:]
                continue
            split_line = line.split()
            pos, __, ref, alts = split_line[1:5]
            position = int(pos) - 1
            if masked:
                if position >= len(mask) or mask[position] == 1:
                    continue
            if interval is not None:
                if position < interval[0] or position >= interval[1]:
                    continue
            split_alts = alts.split(',')
            if "<NON_REF>" in split_alts:
                split_alts.pop(split_alts.index("<NON_REF>"))
            if len(ref) > 1:
                continue  
            if np.any([len(alt) > 1 for alt in split_alts]):
                continue
            if not multiallelic:
                if len(split_alts) > 1:
                    continue
            _gts = [sample.split(':')[0] for sample in split_line[9:]]
            if '.' in "".join(_gts):
                if missing_to_ref:
                    _gts = ['0/0' if '.' in x else x for x in _gts]
                else:
                    warnings.warn("skipping site with missing data")
                    continue
            gts = [re.split("/|\\|", gt) for gt in _gts]           
            _genotypes.append(np.array(gts, dtype=np.int64))
            _sites.append(pos)

    sites = np.array(_sites, dtype=np.int64)
    genotypes = np.array(_genotypes, dtype=np.int64)
    if genotypes.shape == ():
        warnings.warn('Empty genotypes load')
        genotypes = genotypes[:, None]

    return sites, genotypes, sample_ids


def _genotypes_from_str(
    vcf_str, 
    bed_file=None, 
    multiallelic=False,
    missing_to_ref=True,
    interval=None
):
    """

    """
    if bed_file is not None:
        regions = utils.read_bedfile(bed_file)
        mask = utils.regions_to_mask(regions)
        masked = True
    else:
        masked = False

    _sites = []
    _genotypes = []

    for line in vcf_str.split("\n"):
        if line == "":
            continue
        if line.startswith('#'):
            if line.startswith('#CHROM'):
                sample_ids = line.split()[9:]
            continue
        split_line = line.split()
        pos, __, ref, alts = split_line[1:5]
        position = int(pos) - 1
        if masked:
            if position >= len(mask) or mask[position] == 1:
                continue
        if interval is not None:
            if position < interval[0] or position >= interval[1]:
                continue
        split_alts = alts.split(',')
        if "<NON_REF>" in split_alts:
            split_alts.pop(split_alts.index("<NON_REF>"))
        if len(ref) > 1:
            continue  
        if np.any([len(alt) > 1 for alt in split_alts]):
            continue
        if not multiallelic:
            if len(split_alts) > 1:
                continue
        _gts = [sample.split(':')[0] for sample in split_line[9:]]
        if '.' in "".join(_gts):
            if missing_to_ref:
                _gts = ['0/0' if '.' in x else x for x in _gts]
            else:
                warnings.warn("skipping site with missing data")
                continue
        gts = [re.split("/|\\|", gt) for gt in _gts]           
        _genotypes.append(np.array(gts, dtype=np.int64))
        _sites.append(pos)

    sites = np.array(_sites, dtype=np.int64)
    genotypes = np.array(_genotypes, dtype=np.int64)

    return sites, genotypes, sample_ids
