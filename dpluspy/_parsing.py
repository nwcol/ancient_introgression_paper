## Houses functions for computing H2 statistics from genetic data.

from collections import defaultdict
import copy
import gzip
import numpy as np
import re
import scipy
import warnings

from . import utils


"""
## functions for computing H2 statistics from genetic data.
from collections import defaultdict
import copy
import demes
import gzip
import numpy as np
import re
import warnings

from . import utils


def parse_statistics():

    return


def parse_denominators():


    return


def compute_statistics():

    return


def compute_denominators():


    return


def get_H_statistics():

    return


def compute_pairwise_Dp():


    return


def compute_average_pairwise_Dp():


    return


def compute_pairwise_Dp_between():


    return


def compute_average_pairwise_Dp_between():


    return


def _two_haplotype_Dplus():
    ## pairwise

    return


def means_across_region_data():


    return


def get_bootstrap_sets():


    return


def bootstrap_data():


    return


def subset_data():


    return

    

from bisect import bisect
import numpy as np
from h2py import parsing, util


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
"""

def compute_mean_H2(
    data, 
    denominators=None, 
    weighted=False, 
    normalize_to=None
):
    """
    :param normalize_to: normalize to this mutation rate if provided.
        does nothing if `weighted` is not True
    """
    # if `data` is a list, treat its elements as replicate groups of regions
    if isinstance(data, list):
        if denominators is not None:
            for i, rep in enumerate(data):
                for key in rep:
                    assert key in denominators
                    data[i][key] = exchange_denominators(
                        rep[key], denominators[key]
                    )

        means = means_across_replicates(
            data, weighted=weighted, normalize_to=normalize_to
        )

        # compute covariances across replicates
        if len(data) > 1:
            num_reps = len(data)
            num_bins, num_stats = data[0][next(iter(data[0]))]["sums"].shape
            rep_means = np.zeros(
                (num_reps, num_bins, num_stats), dtype=np.float64
            )
            for ii, replicate in enumerate(data):
                rep_means[ii] = means_across_regions(
                    replicate, weighted=weighted, normalize_to=normalize_to
                )
            varcovs = np.zeros(
                (num_bins, num_stats, num_stats), dtype=np.float64
            )
            for ii in range(num_bins):
                varcovs[ii, :, :] = np.cov(rep_means[:, ii, :], rowvar=False)
        else:
            varcovs = None

        example = data[0][next(iter(data[0]))]


    # if `data` is a dict, treat its values as regions
    elif isinstance(data, dict):
        if denominators is not None:
            for key in data:
                data[key] = exchange_denominators(data[key], denominators[key]) 
                
        means = means_across_regions(
            data, weighted=weighted, normalize_to=normalize_to
        ) 

        # compute covariances across regions
        if len(data) > 1:
            num_regions = len(data)
            num_bins, num_stats = data[next(iter(data))]["sums"].shape
            region_means = np.zeros(
                (num_regions, num_bins, num_stats), dtype=np.float64
            )
            for ii, key in enumerate(data):
                region = data[key]
                region_means[ii] = means_across_regions(
                    {'0': region}, weighted=weighted, normalize_to=normalize_to
                )
            varcovs = np.zeros(
                (num_bins, num_stats, num_stats), dtype=np.float64
            )
            for ii in range(num_bins):
                varcovs[ii, :, :] = np.cov(region_means[:, ii, :], rowvar=False)
        else:
            varcovs = None

        example = data[next(iter(data))]

    else:
        raise ValueError("`data` must be list or dict")

    stats = {
        'means': means,
        'covs': varcovs,
        'pops': example['pops'],
        'bins': example['bins']
    }

    return stats


def get_bootstrap_replicates(
    regions, 
    num_reps=None, 
    num_samples=None,
    weighted=False
):
    """

    """
    num_regions = len(regions)

    if num_reps is None: 
        num_reps = num_regions

    if num_samples is None: 
        num_samples = num_regions

    labels = list(regions.keys())
    num_bins, num_stats = regions[labels[0]]["sums"].shape
    replicates = np.zeros((num_reps, num_bins, num_stats), dtype=np.float64)

    for ii in range(num_reps):
        samples = np.random.choice(labels, num_samples, replace=True)
        sampled_regions = {sample: regions[sample] for sample in samples}
        means = means_across_regions(sampled_regions, weighted=weighted)
        replicates[ii] = means

    return replicates


def means_across_regions(regions, weighted=False, normalize_to=None):
    """
    
    """
    num_bins, num_stats = regions[next(iter(regions))]['sums'].shape

    # compute the numerator
    sums = np.zeros((num_bins, num_stats), dtype=np.float64)

    for key in regions:
        sums += regions[key]['sums']

    # compute the denominator with appropriate weighting
    denoms = np.zeros((num_bins), dtype=np.float64)

    if weighted:
        for key in regions:
            if "weights" not in regions[key]:
                raise ValueError("all regions must have precomputed `weights`")

        if normalize_to is None:
            sum_mut = 0  # sum of mean(u) * num_sites across regions
            sum_sites = 0

            for key in regions:
                region = regions[key]
                num_sites = region["weights"]["num_sites"]
                mean_mut = region["weights"]["mean_mut"]
                sum_mut += num_sites * mean_mut
                sum_sites += num_sites

            mean_mut = sum_mut / sum_sites
            mean_mut_sqr = mean_mut ** 2
        else:
            mean_mut_sqr = normalize_to ** 2

        for key in regions:
            mut_prods = regions[key]["weights"]["mut_prods"]
            denoms[:-1] += mut_prods / mean_mut_sqr
            denoms[-1] += regions[key]["denoms"][-1]

    else: 
        for key in regions:
            denoms += regions[key]["denoms"]

    denoms = denoms.reshape((num_bins, 1))
    means = sums / denoms

    if np.any(np.isnan(means)):
        means = np.ma.array(means, mask=np.isnan(means))
    
    return means


def means_across_replicates(replicates, weighted=False, normalize_to=None):
    """
    
    """
    num_reps = len(replicates)
    num_bins, num_stats = replicates[0][next(iter(replicates[0]))]["sums"].shape

    rep_means = np.zeros((num_reps, num_bins, num_stats), dtype=np.float64)

    for i, replicate in enumerate(replicates):
        rep_means[i] = means_across_regions(
            replicate, weighted=weighted, normalize_to=normalize_to
        )

    means = rep_means.mean(0)

    return means



def exchange_denominators(data, denom_data):
    """
    Return a copy of a `data` dictionary after replacing its `denoms` field
    with that of the dictionary `denom_data`.
    """
    ret = copy.deepcopy(data)

    for key in ["denoms", "weights"]:
        if key in denom_data:
            ret[key] = denom_data[key]

    return ret


def enumerate_labels(pops, two_pop=True):
    """
    
    """
    num_pops = len(pops)
    labels = []

    if two_pop:
        for i in range(num_pops):
            for j in range(i, num_pops):
                if i == j:
                    labels.append(tuple(pops[i]))
                else:
                    labels.append((pops[i], pops[j]))
    else:
        for i in range(num_pops):
            labels.append(tuple(pops[i]))

    return labels




def bootstrap_H2(
    regions,
    num_reps=None,
    num_samples=None,
    denom_data=None,
    weighted=False,
    normalize_to=None
):
    """
    
    """
    if num_reps is None: 
        num_reps = len(regions)

    if num_samples is None: 
        num_samples = len(regions)

    keys = list(regions.keys())

    # check that `pops` and `bins` fields match across regions
    for key in keys:
        for field in ['pops', 'bins']:
            if not np.all(regions[key][field] == regions[keys[0]][field]):
                raise ValueError("regions have mismatched `pops` or `bins`")
    
    # replace `denoms` with precomputed denominators if provided
    if denom_data is not None:
        for key in keys:
            assert key in denom_data
            regions[key] = exchange_denominators(regions[key], denom_data[key])
    
    means = means_across_regions(
        regions, weighted=weighted, normalize_to=normalize_to
    )

    num_bins, num_stats = regions[keys[0]]["sums"].shape
    rep_means = np.zeros((num_reps, num_bins, num_stats), dtype=np.float64)

    for ii in range(num_reps):
        samples = np.random.choice(keys, num_samples, replace=True)
        sampled_regions = {sample: regions[sample] for sample in samples}
        rep_mean = means_across_regions(
            sampled_regions, weighted=weighted, normalize_to=normalize_to
        )
        rep_means[ii] = rep_mean

    varcovs = np.zeros((num_bins, num_stats, num_stats), dtype=np.float64)

    for ii in range(num_bins):
        varcovs[ii, :, :] = np.cov(rep_means[:, ii, :], rowvar=False)
    
    data = {
        'pops': regions[keys[0]]['pops'],
        'bins': regions[keys[0]]['bins'],
        'means': means,
        'covs': varcovs,
    }

    return data


## parsing functions


def load_recombination_map(file, map_col="Map(cM)", degree=1):
    """
    Load a recombination map as
    """
    x, y = utils.read_hapmap_rec_map(file, map_col=map_col)
    mapfunc = scipy.interpolate.make_interp_spline(x, y, k=degree)

    return mapfunc


def read_genotypes(
    vcf_file, 
    bed_file=None, 
    region=None,
    min_reg_len=0,
    multiallelic=False,
    missing_to_ref=True,
    report=True,
    return_dict=False
):
    """
    Parse an array of genotypes from a .vcf file. 


    """
    if bed_file is not None:
        regions = utils.read_bedfile(bed_file)
        if min_reg_len:
            regions = regions[(regions[:, 1] - regions[:, 0]) > min_reg_len]
        mask = utils.regions_to_mask(regions)
        masked = True
    else:
        masked = False

    open_func = gzip.open if vcf_file.endswith('.gz') else open

    num_multi = 0
    num_mnv = 0
    num_masked = 0
    data = {}

    with open_func(vcf_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    sample_ids = line.split()[9:]
                continue
            split_line = line.split()
            pos, __, ref, alts = split_line[1:5]
            position = int(pos) - 1
            split_alts = alts.split(',')
            if region is not None:
                if position < region[0] or position >= region[-1]:
                    continue
            # skip the site if it falls outside the mask
            if masked:
                if position >= len(mask) or mask[position] == 1:
                    num_masked += 1
                    continue    
            # this symbol appears in some cases and should be ignored
            if "<NON_REF>" in split_alts:
                split_alts.pop(split_alts.index("<NON_REF>"))
            # check whether the site is a SNP; skip if not
            if len(ref) > 1 or np.any([len(alt) > 1 for alt in split_alts]):
                num_mnv += 1
                continue
            # skip the site if it is multiallelic and `multiallelic` is False
            if not multiallelic:
                if len(split_alts) > 1:
                    num_multi += 1
                    continue
                    
            gt_strs = [sample.split(':')[0] for sample in split_line[9:]]
            # check for missing data
            if '.' in "".join(gt_strs):
                if missing_to_ref:
                    gt_strs = ['0/0' if '.' in x else x for x in gts]
                else:
                    warnings.warn("skipping site with missing data")
                    continue
            gts = [re.split("/|\\|", gt) for gt in gt_strs]           
            line_genotypes = np.array(gts, dtype=np.int64)
            data[position] = line_genotypes

    if return_dict:
        ret = (data, sample_ids)
    else:
        sites = np.array(list(data.keys()), dtype=np.int64)
        genotypes = np.array([data[pos] for pos in data], dtype=np.int64)
        ret = (sites, genotypes, sample_ids)
    if report:
        pass

    return ret


def parse_statistics(
    vcf_file,
    region=None,
    rec_map_file=None,
    uniform_r=None,
    mut_map_file=None,
    r_bins=None,
    bp_bins=None,
    bed_file=None,
    pop_file=None,
    use_haplotypes=False,
    compute_two_pop=True,
    compute_denom=True,
    snp_denom=False,
    min_reg_len=0,
    L=None,
    map_degree=1
):
    """
    
    """
    sites, genotypes, sample_ids = read_genotypes(
        vcf_file, bed_file=bed_file, min_reg_len=min_reg_len
    )

    if rec_map_file is not None:
        site_map = utils.get_rec_map(rec_map_file, sites)
        use_r_bins = True
    elif uniform_r is not None:
        site_map = utils.get_uniform_rec_map(uniform_r, sites)
        use_r_bins = True
    else:
        site_map = sites
        use_r_bins = False
        warnings.warn('binning site pairs by physical distance')
    
    if use_r_bins:
        if r_bins is not None:
            if np.any(np.diff(r_bins)) <= 0:
                raise ValueError("bins are not monotonic")
            ret_bins = r_bins
        else:
            raise ValueError("you must provide bins")
        bins = utils.map_function(ret_bins)
    else:
        if bp_bins is not None:
            bins = bp_bins
        else:
            raise ValueError("you must provide bins")
        ret_bins = bins

    stats = compute_statistics(
        genotypes,
        sample_ids,
        site_map,
        bins,
        ret_bins=ret_bins,
        region=region,
        sites=sites,
        pop_file=pop_file,
        use_haplotypes=use_haplotypes,
        compute_two_pop=compute_two_pop
    )

    if snp_denom:
        compute_denom = True

    if compute_denom:
        if snp_denom:
            positions = sites

        else:
            if bed_file is None:
                if L is not None:
                    bed_regions = np.array([[0, L]])
                else: 
                    raise ValueError('you must provide `L` or `bed_file`')
            else:
                bed_regions = utils.read_bedfile(bed_file)

            if min_reg_len:
                filt = (bed_regions[:, 1] - bed_regions[:, 0]) > min_reg_len
                bed_regions = bed_regions[filt]

            if region is not None:
                start, end = region[0], region[-1]
            else:
                start, end = bed_regions[0, 0], bed_regions[-1, 1]

            # mask = util.regions_to_mask(bed_regions)[start:end]
            # positions = np.nonzero(mask)[0] + start

            positions = np.nonzero(~utils.regions_to_mask(bed_regions))[0]
            positions = positions[(positions >= start) & (positions < end)]

        if rec_map_file is not None:
            pos_map = utils.get_rec_map(rec_map_file, positions)
        elif uniform_r is not None:
            pos_map = utils.get_uniform_rec_map(uniform_r, positions)
        else:
            pos_map = positions

        if mut_map_file is not None:
            mut_map = utils.read_mutation_map(mut_map_file, positions)
            stats["denoms"], stats["weights"] = compute_denominators(
                pos_map, 
                bins, 
                region=region, 
                positions=positions,
                mut_map=mut_map
            )
        else:
            stats["denoms"] = compute_denominators(
                pos_map, bins, region=region, positions=positions
            )

    return stats


def compute_statistics(   
    genotypes,
    sample_ids,
    site_map,
    bins,
    ret_bins=None,
    region=None,
    sites=None,
    pop_file=None,
    pop_dict=None,
    use_haplotypes=False,
    compute_two_pop=True
):
    """
    
    """
    bins = np.asanyarray(bins)

    if region is not None:
        if sites is None:
            raise ValueError('you must provide sites to subset to a region')
        if region.shape == (3,):
            start, right_lim = np.searchsorted(sites, [region[0], region[2]])
            left_lim = np.searchsorted(sites[start:], region[1])
        elif region.shape == (2,):
            start, right_lim = np.searchsorted(sites, region)
            left_lim = len(site_map)
        else:
            raise ValueError('invalid region shape')
        genotypes = genotypes[start:right_lim]
        site_map = site_map[start:right_lim]
    else:
        left_lim = len(site_map)

    if pop_file is not None: 
        pop_dict = defaultdict(list)
        with open(pop_file, 'r') as popf:
            for line in popf:
                sample, pop = line.split()
                pop_dict[pop].append(sample)

    if pop_dict is not None:
        pops = list(pop_dict.keys())
        pop_indices = {}
        for pop in pops:
            samples = pop_dict[pop]
            indices = [sample_ids.index(sample) for sample in samples]
            pop_indices[pop] = np.array(indices)
    else:
        pops = sample_ids
        pop_indices = {pop: [i] for i, pop in enumerate(pops)}

    pop_genotypes = {pop: genotypes[:, pop_indices[pop]] for pop in pops}

    short_genotypes = {pop: pop_genotypes[pop][:left_lim] for pop in pops}
    sums_H = compute_h_stat_sums(
        short_genotypes, 
        compute_two_pop=compute_two_pop
    )
    sums_H2 = compute_h2_stat_sums(
        pop_genotypes,
        site_map,
        bins,
        compute_two_pop=compute_two_pop,
        use_haplotypes=use_haplotypes,
        left_lim=left_lim
    )
    sums = np.vstack((sums_H2, sums_H[np.newaxis, :]))

    stats = {
        "bins": ret_bins if ret_bins is not None else bins,
        "pops": pops,
        "sums": sums
    }

    return stats


def compute_h_stat_sums(pop_genotypes, compute_two_pop=True, verbose=True):

    pops = list(pop_genotypes.keys())
    num_pops = len(pops)

    if compute_two_pop:
        num_stats = num_pops + num_pops * (num_pops - 1) // 2
        pop_indices = [(i, j) for i in range(num_pops) 
                       for j in range(i, num_pops)]
    else:
        num_stats = num_pops
        pop_indices = [(i, i) for i in range(num_pops)]

    H_sums = np.zeros(num_stats, dtype=np.float64)

    for k, (i, j) in enumerate(pop_indices):
        if i == j:
            pop = pops[i]
            H_sums[k] = _one_pop_h(pop_genotypes[pop])
        else:
            pop1 = pops[i]
            pop2 = pops[j]
            H_sums[k] = _two_pop_h(pop_genotypes[pop1], pop_genotypes[pop2])

    if verbose:
        num_sites = len(pop_genotypes[next(iter(pop_genotypes))])
        print(utils.get_time(), f"computed H sums for {num_sites} sites")

    return H_sums


def _one_pop_h(genotypes):

    num_sites, _n, _ = genotypes.shape
    n = 2 * _n
    variants = np.reshape(genotypes, (num_sites, n))

    numer = 0
    for i in range(n - 1):
        for j in range(i + 1, n):
            numer += (variants[:, i] != variants[:, j]).sum()

    denom = n * (n - 1) / 2
    sum_H = numer / denom
    return sum_H


def _two_pop_h(genotypes1, genotypes2):

    num_sites1, _n1, _ = genotypes1.shape
    num_sites2, _n2, _ = genotypes2.shape
    n1, n2 = 2 * _n1,  2 * _n2
    variants1 = np.reshape(genotypes1, (num_sites1, n1))
    variants2 = np.reshape(genotypes2, (num_sites2, n2))

    numer = 0
    for i in range(n1):
        for j in range(n2):
            numer += (variants1[:, i] != variants2[:, j]).sum()

    denom = n1 * n2
    sum_H = numer / denom
    return sum_H


def compute_h2_stat_sums(
    pop_genotypes,
    site_map,
    bins,
    compute_two_pop=True,
    use_haplotypes=False,
    left_lim=None,
    verbose=True
):
    """
    
    """
    pops = list(pop_genotypes.keys())
    
    if use_haplotypes:
        pop_haplotypes = {}
        for pop in pop_genotypes:
            genotypes = pop_genotypes[pop]
            num_sites, n, _ = genotypes.shape
            pop_haplotypes[pop] = np.reshape(genotypes, (num_sites, 2 * n))

    num_bins = len(bins) - 1
    num_pops = len(pops)

    if compute_two_pop:
        num_stats = num_pops + num_pops * (num_pops - 1) // 2
        pop_indices = [(i, j) for i in range(num_pops) 
                       for j in range(i, num_pops)]
    else:
        num_stats = num_pops
        pop_indices = [(i, i) for i in range(num_pops)]

    H2_sums = np.zeros((num_bins, num_stats), dtype=np.float64)

    for k, (i, j) in enumerate(pop_indices):
        if i == j:
            pop = pops[i]

            if use_haplotypes:
                H2_sums[:, k] = haplotype_h2(
                    pop_haplotypes[pop],
                    site_map,
                    bins,
                    left_lim=left_lim
                )
            else:
                H2_sums[:, k] = genotype_h2(
                    pop_genotypes[pop],
                    site_map,
                    bins,
                    left_lim=left_lim
                )
        else:
            pop1 = pops[i]
            pop2 = pops[j]

            if use_haplotypes:
                H2_sums[:, k] = two_pop_haplotype_h2(
                    pop_haplotypes[pop1],
                    pop_haplotypes[pop2],
                    site_map,
                    bins,
                    left_lim=left_lim
                )
            else:
                H2_sums[:, k] = two_pop_genotype_h2(
                    pop_genotypes[pop1],
                    pop_genotypes[pop2],
                    site_map,
                    bins,
                    left_lim=left_lim
                )
    
    if verbose:
        print(utils.get_time(), f"computed H2 sums for {left_lim} left loci")

    return H2_sums


## functions for computing single/multi-sample H2


def haplotype_h2(
    haplotypes, 
    site_map, 
    bins, 
    left_lim=None
):
    """
    Compute one-population haplotype H2 by averaging binned H2 sums across
    (2 * n) choose 2 within-population haplotype pairs.
    """
    _, n = haplotypes.shape
    num_bins = len(bins) - 1

    if n == 2:
        num_H2 = _two_haplotype_h2(
            haplotypes[:, 0],
            haplotypes[:, 1],
            site_map,
            bins, 
            left_lim=left_lim
        )
    else:
        pair_H2 = np.zeros((num_bins, n * (n - 1) // 2), dtype=np.float64)
        indices = [(i, j) for i in range(n) for j in range(i + 1, n)]

        for k, (i, j) in enumerate(indices):
                pair_H2[:, k] = _two_haplotype_h2(
                    haplotypes[:, i],
                    haplotypes[:, j],
                    site_map,
                    bins, 
                    left_lim=left_lim
                )
        num_H2 = pair_H2.mean(1)

    return num_H2


def two_pop_haplotype_h2(
    haplotypes1,
    haplotypes2,
    site_map,
    bins,
    left_lim=None
):
    """
    Compute two-population haplotype H2 as the mean of H2 across 2 * n1 * n2
    cross-population haplotype pairs. Returns binned average sums.
    """
    num_sites1, n1 = haplotypes1.shape
    num_sites2, n2 = haplotypes2.shape
    num_bins = len(bins) - 1

    assert num_sites1 == num_sites2

    pair_H2 = np.zeros((num_bins, n1 * n2), dtype=np.float64)
    indices = [(i, j) for i in range(n1) for j in range(n2)]

    for k, (i, j) in enumerate(indices):
            pair_H2[:, k] = _two_haplotype_h2(
                haplotypes1[:, i],
                haplotypes2[:, j],
                site_map,
                bins, 
                left_lim=left_lim
            )
    num_H2 = pair_H2.mean(1)

    return num_H2


def genotype_h2(
    genotypes,
    site_map,
    bins,
    left_lim=None
):
    """
    Compute one-population genotype H2. If there is more than one sample in the
    population, then we take the average of single-sample genotype H2.
    """
    _, n, __ = genotypes.shape
    num_bins = len(bins) - 1

    if n == 1:
        num_H2 = _one_genotype_h2(
            genotypes[:, 0],
            site_map,
            bins, 
            left_lim=left_lim
        )
    else:
        sample_H2 = np.zeros((num_bins, n), dtype=np.float64)

        for i in range(n):
                sample_H2[:, i] = _one_genotype_h2(
                    genotypes[:, i],
                    site_map,
                    bins, 
                    left_lim=left_lim
                )
        num_H2 = sample_H2.mean(1)

    return num_H2


def two_pop_genotype_h2(
    genotypes1,
    genotypes2,
    site_map,
    bins,
    left_lim=None    
):
    """
    Compute two-population genotype H2 by taking the mean across the n1 * n2
    cross-population sample pairs. Returns the statistic as binned sums.
    """
    num_sites1, n1, _ = genotypes1.shape
    num_sites2, n2, _ = genotypes2.shape
    num_bins = len(bins) - 1

    assert num_sites1 == num_sites2

    pair_H2 = np.zeros((num_bins, n1 * n2), dtype=np.float64)
    indices = [(i, j) for i in range(n1) for j in range(n2)]

    for k, (i, j) in enumerate(indices):
            pair_H2[:, k] = _two_genotype_h2(
                genotypes1[:, i],
                genotypes2[:, j],
                site_map,
                bins, 
                left_lim=left_lim
            )
    num_H2 = pair_H2.mean(1)
    return num_H2


## lowest-level functions to compute H2 for single samples/sample pairs


def _one_genotype_h2(
    genotypes,
    site_map,
    bins,
    left_lim=None
):
    """
    A wrapper for `_two_haplotype_h2`. Computes binned sums of H2 for a single
    array of genotypes. 
    """
    num_H2 = _two_haplotype_h2(
        genotypes[:, 0],
        genotypes[:, 1],
        site_map,
        bins,
        left_lim=left_lim
    )

    return num_H2


def _two_haplotype_h2(
    haplotype1, 
    haplotype2,
    site_map,  
    bins,  
    left_lim=None
):
    """
    Compute bin-sum H2 between two arrays of haplotypes.
    """
    assert len(haplotype1) == len(haplotype2)
    assert len(haplotype1) == len(site_map)

    indicator = haplotype1 != haplotype2
    num_H2 = _bin_pair_products(indicator, site_map, bins, left_lim=left_lim)

    return num_H2


def _two_genotype_h2(
    genotypes1,
    genotypes2,
    site_map,
    bins,
    left_lim=None
):  
    """
    Compute bin-sum H2 between two arrays of one-sample genotypes.
    """
    assert len(genotypes1) == len(genotypes2)
    assert len(genotypes1) == len(site_map)   
 
    indicator = genotypes1[:, :, np.newaxis] != genotypes2[:, np.newaxis]
    site_weights = indicator.sum((2, 1)) / 4
    num_H2 = _bin_pair_products(site_weights, site_map, bins, left_lim=left_lim)
    
    return num_H2


## function for computing denominators of H2, H statistics


def compute_denominators(
    pos_map, 
    bins, 
    positions=None,
    region=None,
    mut_map=None,
    verbose=True
):
    """
    
    """
    if region is not None:
        if positions is None:
            raise ValueError('you must provide positions to subset to a region')
        if region.shape == (3,):
            start, right_lim = np.searchsorted(positions, [region[0],region[2]])
            left_lim = np.searchsorted(positions[start:], region[1])
        elif region.shape == (2,):
            start, right_lim = np.searchsorted(positions, region)
            left_lim = None
        else:
            raise ValueError('invalid region shape')

        pos_map = pos_map[start:right_lim]

        if mut_map is not None:
            mut_map = mut_map[start:right_lim]
    else:
        left_lim = len(pos_map)

    denom_H = left_lim
    denoms_H2 = _fast_bin_pair_counts(pos_map, bins, left_lim=left_lim)
    denoms = np.append(denoms_H2, denom_H)

    if mut_map is None:
        ret = denoms
    else:
        mut_prods = _fast_bin_pair_products(
            mut_map, pos_map, bins, left_lim=left_lim
        )
        weights = {
            "mut_prods": mut_prods,
            "num_sites": int(left_lim),
            "mean_mut": float(mut_map.mean())
        }
        ret = (denoms, weights)

    if verbose:
        print(utils.get_time(), f"computed denominator for {left_lim} left loci")

    return ret


## functions for computing locus pair sums


def _bin_pair_counts(site_map, bins, left_lim=None):
    """
    
    """
    if not left_lim:
        left_lim = len(site_map)

    verbose = 1e6

    cum_nums = np.zeros(len(bins), dtype=int)

    for i, rl in enumerate(site_map[:left_lim]):
        edges = np.searchsorted(site_map[i + 1:], rl + bins)
        cum_nums += edges

        if i % verbose == 0 and i > 0:
            print(utils.get_time(), f'num pairs counted at site {i}')

    num_pairs = np.diff(cum_nums)

    return num_pairs


def _fast_bin_pair_counts(rec_map, bins, left_lim=None, verbose=True):
    """
    
    """
    if not left_lim:
        left_lim = len(rec_map)

    num_bins = len(bins) - 1
    num_pairs = np.zeros(num_bins, dtype=int)

    if bins[0] == 0:
        edge0 = np.arange(1, left_lim + 1)
    else:
        edge0 = np.searchsorted(rec_map, rec_map[:left_lim] + bins[0])

    for i, b in zip(range(num_bins), bins[1:]):
        edge1 = np.searchsorted(rec_map, rec_map[:left_lim] + b)
        num_pairs[i] = (edge1 - edge0).sum() 
        edge0 = edge1
        
        if verbose:
            print(utils.get_time(), f"site pair counts computed in bin {i}")

    return num_pairs


def _bin_pair_products(site_vals, site_map, bins, left_lim=None, verbose=1e6):
    # get binned sums of site products
    if len(site_vals) != len(site_map):
        raise ValueError("map/value length mismatch")
    
    if left_lim is None:
        left_lim = len(site_map)

    if verbose is None:
        verbose = 1e10

    cum_vals = np.cumsum(site_vals)
    cum_prods = np.zeros(len(bins), dtype=np.float64)

    for i, (rl, ul) in enumerate(zip(site_map[:left_lim], site_vals[:left_lim])):
        if ul > 0:
            edges = np.searchsorted(site_map[i + 1:], rl + bins)
            cum_prods += ul * cum_vals[i:][edges]

            if i % verbose == 0 and i > 0:
                print(utils.get_time(), f"site pair products computed at locus {i}")

    prod_sums = np.diff(cum_prods)

    return prod_sums


def _fast_bin_pair_products(
    site_vals, 
    site_map, 
    bins, 
    left_lim=None,
    verbose=True
):
    """
    the bin-vectorized version (fast but large memory expense)
    """
    if not left_lim:
        left_lim = len(site_map)

    cum_vals = np.cumsum(site_vals)
    left_vals = site_vals[:left_lim]

    num_bins = len(bins) - 1
    products = np.zeros(num_bins, dtype=float)

    if bins[0] == 0:
        indices = np.arange(0, left_lim)
    else:
        indices = np.searchsorted(site_map, site_map[:left_lim] + bins[0]) - 1

    cum_sum0 = cum_vals[indices]

    for i, b in zip(range(num_bins), bins[1:]):
        indices = np.searchsorted(site_map, site_map[:left_lim] + b) - 1
        cum_sum1 = cum_vals[indices]
        cum_products = left_vals * (cum_sum1 - cum_sum0)
        products[i] = cum_products.sum()
        cum_sum0 = cum_sum1

        if verbose:
            print(utils.get_time(), f"site pair products computed in bin {i}")

    return products


## sampling functions for numerically computing H2


def _sample_haplotype_h2(haplotypes, num_reps=1000):
    # from an array of haplotypes with shape (2, n), sample haplotypes to 
    # compute expected H2 numerically. for numerical validation
    num_sites, num_haps = haplotypes.shape
    assert num_sites == 2

    indices = np.arange(num_haps)
    sum_H2 = 0
    for i in range(num_reps):
        sample = haplotypes[:, np.random.choice(indices, size=2, replace=False)]
        sum_H2 += (
            (sample[0, 0] != sample[0, 1]) 
            & (sample[1, 0] != sample[1, 1])
        )
    
    mean_H2 = sum_H2 / num_reps
    return mean_H2


def _sample_two_pop_haplotype_h2(haplotypes1, haplotypes2, num_reps=1000):
    # for numerical validation
    num_sites1, num_haps1 = haplotypes1.shape
    num_sites2, num_haps2 = haplotypes2.shape
    assert num_sites1 == num_sites2 == 2

    sum_H2 = 0
    for i in range(num_reps):
        sample1 = haplotypes1[:, np.random.randint(num_haps1)]
        sample2 = haplotypes2[:, np.random.randint(num_haps2)]
        sum_H2 += ((sample1[0] != sample2[0]) & (sample1[1] != sample2[1]))
    
    mean_H2 = sum_H2 / num_reps
    return mean_H2


def _sample_genotype_h2(genotypes, num_reps=1000, between=False):
    #
    num_sites, num_samps, _ = genotypes.shape
    assert num_sites == 2

    # precompute H2 for each sample
    sample_H2s = []
    for i in range(num_samps):
        genotype = genotypes[:, i]
        sample_H2s.append(
            (genotype[0, 0] != genotype[0, 1]) 
            & (genotype[1, 0] != genotype[1, 1]) 
        )
    sample_H2s = np.array(sample_H2s)

    sum_H2 = 0

    # allow haplotypes to be sampled from different genomes
    if between:
        for i in range(num_reps):
            index1 = np.random.randint(num_samps)
            index2 = np.random.randint(num_samps)

            if index1 == index2:
                sum_H2 += sample_H2s[index1]
            else:
                sample1 = genotypes[:, index1]
                sample2 = genotypes[:, index2]
                hap1 = sample1[[0, 1], np.random.randint(2, size=2)]
                hap2 = sample2[[0, 1], np.random.randint(2, size=2)]   
                sum_H2 += ((hap1[0] != hap2[0]) & (hap1[1] != hap2[1]))    

    # average over within-sample H2
    else:
        for i in range(num_reps):
            sum_H2 += sample_H2s[np.random.randint(num_samps)]

    mean_H2 = sum_H2 / num_reps
    return mean_H2


def _sample_two_pop_genotype_h2(genotypes1, genotypes2, num_reps=1000):
    # for numerical validation. two-sample genotype H2 is the simpler case!
    num_sites1, num_samps1, _ = genotypes1.shape
    num_sites2, num_samps2, _ = genotypes2.shape
    assert num_sites1 == num_sites2 == 2

    sum_H2 = 0
    for i in range(num_reps):
        sample1 = genotypes1[:, np.random.randint(num_samps1)]
        sample2 = genotypes2[:, np.random.randint(num_samps2)]
        hap1 = sample1[[0, 1], np.random.randint(2, size=2)]
        hap2 = sample2[[0, 1], np.random.randint(2, size=2)]
        sum_H2 += ((hap1[0] != hap2[0]) & (hap1[1] != hap2[1]))
    
    mean_H2 = sum_H2 / num_reps
    return mean_H2


## probably obsolete but useful for reference


def tally_haplotypes(h1, h2):
    #
    n = len(h1)
    n11 = (h1 & h2).sum()
    n10 = h1.sum() - n11
    n01 = h2.sum() - n11
    n00 = n - n11 - n10 - n01
    return (n11, n10, n01, n00)


def tally_genotypes(g1, g2):
    # tallies up two-locus genotypes
    # note that & is a symbol for np.logical_and()
    n = len(g1)
    n22 = ((g1 == 2) & (g2 == 2)).sum()
    n21 = ((g1 == 2) & (g2 == 1)).sum() 
    n20 = (g1 == 2).sum() - n22 - n21
    n12 = ((g1 == 1) & (g2 == 2)).sum()
    n11 = ((g1 == 1) & (g2 == 1)).sum()
    n10 = (g1 == 1).sum() - n12 - n11
    n02 = ((g1 == 0) & (g2 == 2)).sum()
    n01 = ((g1 == 0) & (g2 == 1)).sum()
    n00 = n - n22 - n21 - n20 - n12 - n11 - n10 - n02 - n01
    return (n22, n21, n20, n12, n11, n10, n02, n01, n00)


def _haplotype_h2_from_counts(counts):
    #
    c1, c2, c3, c4 = counts
    numer = c1 * c4 + c2 * c3
    num = c1 + c2 + c3 + c4
    h2 = numer / (num * (num - 1) / 2)    
    return h2


def _two_pop_haplotype_h2_from_counts(counts1, counts2):
    #
    c11, c12, c13, c14 = counts1 
    c21, c22, c23, c24 = counts2 
    numer = c11 * c24 + c14 * c21 + c12 * c23 + c13 * c22
    num1 = c11 + c12 + c13 + c14
    num2 = c21 + c22 + c23 + c24
    h2 = numer / (num1 * num2)
    return h2


def _genotype_h2_from_counts(counts):
    # shape (9, b). counts of two-locus genotypes 1/1-1/1, ... 0/0-0/0
    n1, n2, n3, n4, n5, n6, n7, n8, n9 = counts
    numer = (
        n1 * n5 / 4
        + n1 * n6 / 2
        + n1 * n8 / 2
        + n1 * n9
        + n2 * n4 / 4
        + n2 * n5 / 4
        + n2 * n6 / 4
        + n2 * n7 / 2
        + n2 * n8 / 4
        + n2 * n9 / 2
        + n3 * n4 / 2
        + n3 * n5 / 4
        + n3 * n7
        + n3 * n8 / 2
        + n4 * n5 / 4
        + n4 * n6 / 2
        + n4 * n8 / 4
        + n4 * n9 / 2
        + n5  ######
        + n5 * n6 / 4 
        + n5 * n7 / 4
        + n5 * n8 / 4
        + n5 * n9 / 4
        + n6 * n7 / 2
        + n6 * n8 / 4
    )
    num = counts.sum(0) 
    h2 = numer / (num * (num - 1) / 2)  
    return h2


def _two_pop_genotype_h2_from_counts(counts1, counts2):
    # shapes (9, b)
    n11, n12, n13, n14, n15, n16, n17, n18, n19 = counts1
    n21, n22, n23, n24, n25, n26, n27, n28, n29 = counts2
    numer = (
        n11 * n29 + n19 * n21
        + (n11 * n26 + n16 * n21) / 2
        + (n11 * n28 + n18 * n21) / 2
        + (n11 * n29 + n19 * n21) / 2
        + (n12 * n24 + n14 * n22) / 4
        + (n12 * n25 + n15 * n22) / 4
        + (n12 * n26 + n16 * n22) / 4
        + (n12 * n27 + n17 * n22) / 2
        + (n12 * n28 + n18 * n22) / 4
        + (n12 * n29 + n19 * n22) / 2
        + (n13 * n24 + n14 * n23) / 2
        + (n13 * n25 + n15 * n23) / 4
        + n13 * n27 + n17 * n23
        + (n13 * n28 + n18 * n23) / 2
        + (n14 * n25 + n15 * n24) / 4
        + (n14 * n26 + n16 * n24) / 2
        + (n14 * n28 + n18 * n24) / 4
        + (n14 * n29 + n19 * n24) / 2
        + (n15 * n25) / 4 
        + (n15 * n26 + n16 * n25) / 4 
        + (n15 * n27 + n17 * n25) / 4
        + (n15 * n28 + n19 * n25) / 4
        + (n15 * n29 + n19 * n25) / 4
        + (n16 * n27 + n17 * n26) / 2
        + (n16 * n28 + n18 * n26) / 4
    )
    num1 = counts1.sum(0)
    num2 = counts1.sum(0)
    h2 = numer / (num1 * num2)
    return h2
