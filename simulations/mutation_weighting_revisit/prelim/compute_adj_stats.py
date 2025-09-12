
import dpluspy 
import numpy as np
import pickle
import sys


def main():

    vcf_fname = sys.argv[1]
    bed_fname = sys.argv[2]
    rec_fname = sys.argv[3]
    mut_fname = sys.argv[4]
    interval_fname = sys.argv[5]
    chrom = int(sys.argv[6])
    u_bar = float(sys.argv[7])
    out_fname = sys.argv[8]
   
    r_bins = np.logspace(-6, -2, 17)

    adj_sums = parse_stats(
        vcf_fname,
        bed_file=bed_fname,
        rec_map_file=rec_fname,
        mut_map_file=mut_fname,
        interval_file=interval_fname,
        pos_col="pos",
        map_col="cM",
        mut_col="mut_map",
        r_bins=r_bins,
        chrom=chrom,
        u_bar=u_bar
    )

    with open(out_fname, "wb") as fout:
        pickle.dump(adj_sums, fout)

    return


def parse_stats(
    vcf_file,
    u_bar=1e-8,
    ts_sample_ids=None,
    bed_file=None,
    pop_file=None,
    pop_mapping=None,
    rec_map_file=None,
    pos_col="Position(bp)",
    map_col="Map(cM)",
    map_sep=None,
    interp_method="linear",
    r=None,
    r_bins=None,
    bp_bins=None,
    mut_map_file=None,
    mut_col=None,
    interval=None,
    intervals=None,
    interval_file=None,
    chrom="None",
    phased=False,
    get_cross_pop=True,
    get_denoms=True,
    allow_multi=True,
    missing_to_ref=False,
    apply_filter=False,
    verbose=True
):
    if interval is not None:
        intervals = [interval]
    if interval_file is not None:
        intervals = np.loadtxt(interval_file)
    if np.array(intervals).ndim == 1:
        intervals = np.array(intervals)[None, :]
    # Convert intervals into a list of 1d arrays
    intervals = [np.asarray(x).flatten().astype(np.int64) for x in intervals]

    if get_denoms:
        if bed_file is not None: 
            positions = dpluspy.utils._read_bed_file_positions(bed_file) + 1
            seq_length = positions[-1]
        else:
            seq_length = intervals[-1][-1]
            positions = np.arange(1, seq_length)

        if mut_map_file is not None:
            mut_map = dpluspy.parsing._load_mutation_map(
                mut_map_file, positions, map_col=mut_col)
        else:
            raise ValueError("You need to provide a mutation map!")
    else:
        positions = None
        mut_map = None
        seq_length = intervals[-1][-1] + 1

    if r_bins is not None:
        if isinstance(r_bins, str):
            r_bins = np.loadtxt(r_bins)
        # Convert bins in r to Morgans 
        bins = dpluspy.utils._map_function(r_bins)
        if rec_map_file is not None:
            map_fxn = dpluspy.parsing._load_recombination_map(
                rec_map_file, 
                pos_col=pos_col,
                map_col=map_col,
                interp_method=interp_method,
                map_sep=map_sep
            )
        elif r is not None:
            map_fxn = \
                dpluspy.parsing._get_uniform_recombination_map(r, seq_length)
        else:
            raise ValueError("You must provide recombination map information")
        # Save r bins for output
        ret_bins = r_bins
    elif bp_bins is not None:
        map_fxn = lambda x: x
        ret_bins = None
    else:
        raise ValueError("You must provide bins")

    if pop_file is not None:
        pop_mapping = load_label_file(pop_file)

    if pop_mapping is not None:
        vcf_ids = [sample for pop in pop_mapping for sample in pop_mapping[pop]]
    else:
        vcf_ids = None

    # Read genotypes from a VCF file or extract them from a tree sequence
    if isinstance(vcf_file, str):
        sites, genotypes, sample_ids = dpluspy.parsing.get_vcf_genotypes(
            vcf_file, 
            sample_ids=vcf_ids,
            bed_file=bed_file, 
            allow_multi=allow_multi,
            missing_to_ref=missing_to_ref,
            apply_filter=apply_filter
        )
    else:
        sites, genotypes, sample_ids = dpluspy.parsing.get_ts_genotypes(
            vcf_file, 
            ts_sample_ids=ts_sample_ids,
            sample_ids=vcf_ids,
            bed_file=bed_file, 
            allow_multi=allow_multi,
            missing_to_ref=missing_to_ref,
            apply_filter=apply_filter
        )
    # Construct a dict mapping population IDs to population genotype arrays
    genotype_dict = get_genotype_dict(
        genotypes, sample_ids, sample_labels=pop_mapping)
    
    stats = compute_stats(    
        sites,
        genotype_dict,
        map_fxn,
        bins,
        intervals,
        positions=positions,
        mut_map=mut_map,
        u_bar=u_bar,
        chrom=chrom,
        get_cross_pop=get_cross_pop,
        phased=phased,
        ret_bins=ret_bins,
        verbose=verbose
    )
    return stats


def compute_stats(
    sites,
    genotype_dict,
    map_func,
    bins,
    intervals,
    u_bar=None,
    positions=None,
    mut_map=None,
    chrom="None",
    get_cross_pop=True,
    phased=False,
    verbose=True,
    ret_bins=None
):

    samples = list(genotype_dict.keys())
    ret = dict()

    for ii, interval in enumerate(intervals):
        assert len(interval) == 3
        left_interval = interval[:2]
        right_interval = interval[1:]

        stats = dict()
        stats["bins"] = ret_bins
        stats["pop_ids"] = samples

        if positions is not None:
            stats["denoms"] = dpluspy.parsing.denoms_within(
                positions, map_func, bins, left_interval)
        stats["sums"] = get_stats_within(
            sites, 
            left_interval, 
            genotype_dict, 
            map_func, 
            bins, 
            mut_map,
            u_bar,
            get_cross_pop=get_cross_pop
        )
        if verbose:
            print(dpluspy.utils._current_time(), 
                f"Computed stats within chrom {chrom} interval {ii} "
                f"{interval[0]}:{interval[1]}")
        
        if right_interval[1] > right_interval[0]:
            if positions is not None:
                stats["denoms"] += dpluspy.parsing.denoms_between(
                    positions, map_func, bins, 
                    (left_interval, right_interval))
            stats["sums"] += get_stats_between(
                sites, 
                (left_interval, right_interval), 
                genotype_dict, 
                map_func, 
                bins,
                mut_map,
                u_bar,
                get_cross_pop=get_cross_pop
            )
            if verbose:
                print(dpluspy.utils._current_time(), 
                    f"Computed stats between chrom {chrom} intervals {ii} "
                    f"{left_interval[0]}:{right_interval[1]}:{interval[2]}")

        key = (chrom, ii)
        ret[key] = stats
    
    return ret


def get_stats_within(
    sites, 
    interval,
    genotype_dict, 
    map_func, 
    bins,
    mut_map, 
    u_bar,
    get_cross_pop=True,
):
    """
    """
    start, end = interval
    where = np.where((sites >= start) & (sites < end))[0]
    sub_genotype_dict = {p: genotype_dict[p][where] for p in genotype_dict}
    rec_map = map_func(sites[where])
    mut_map = mut_map[where]
    sums = compute_stats_within(
        sub_genotype_dict, 
        rec_map, 
        bins,
        mut_map,
        u_bar,
        cross_pop=get_cross_pop,
    )
    return sums


def get_stats_between(    
    sites, 
    intervals,
    genotype_dict, 
    map_func, 
    bins, 
    mut_map, 
    u_bar,
    get_cross_pop=True
):
    """
    Higher-level than `compute_stats_within`. Subsets loaded data 
    """
    (left_start, left_end), (right_start, right_end) = intervals
    where_left = np.where((sites >= left_start) & (sites < left_end))[0]
    left_genotype_dict = {
        p: genotype_dict[p][where_left] for p in genotype_dict}
    left_rec_map = map_func(sites[where_left])
    left_mut_map = mut_map[where_left]

    where_right = np.where((sites >= right_start) & (sites < right_end))[0]
    right_genotype_dict = {
        p: genotype_dict[p][where_right] for p in genotype_dict}
    right_rec_map = map_func(sites[where_right])
    right_mut_map = mut_map[where_right]

    sums = compute_stats_between(
        left_genotype_dict,
        right_genotype_dict,
        left_rec_map,
        right_rec_map,
        bins, 
        left_mut_map,
        right_mut_map,
        u_bar,
        cross_pop=get_cross_pop
    )
    return sums


def compute_stats_within(
    genotype_dict, 
    rec_map, 
    bins,
    mut_map,
    u_bar,
    cross_pop=True,
):

    pop_ids = list(genotype_dict.keys())
    num_pops = len(pop_ids)
    if cross_pop:
        num_stats = (num_pops + num_pops ** 2) // 2
    else:
        num_stats = num_pops
    sums = np.zeros((len(bins), num_stats))
    idx = 0
    for ii, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[ii:]:
            if pop_i == pop_j:
                Gt_ii = genotype_dict[pop_i]
                sums[:-1, idx] = compute_numer_one_pop_within(
                    Gt_ii, rec_map, bins, mut_map, u_bar)
            else:
                if not cross_pop:
                    continue
                Gt_ii = genotype_dict[pop_i]
                Gt_jj = genotype_dict[pop_j]
                sums[:-1, idx] = compute_numer_cross_pop_within(
                    Gt_ii, Gt_jj, rec_map, bins, mut_map, u_bar)
            idx += 1
    # sums[-1] = dpluspy.parsing._compute_pi(genotype_dict, cross_pop=cross_pop)
    return sums


def compute_stats_between(
    left_genotype_dict,
    right_genotype_dict,
    left_rec_map,
    right_rec_map,
    bins, 
    left_mut_map,
    right_mut_map,
    u_bar,
    cross_pop=True
):

    pop_ids = list(left_genotype_dict.keys())
    num_pops = len(pop_ids)
    if cross_pop:
        num_stats = (num_pops + num_pops ** 2) // 2
    else:
        num_stats = num_pops
    sums = np.zeros((len(bins), num_stats))
    idx = 0
    for ii, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[ii:]:
            if pop_i == pop_j:
                left_Gt_i = left_genotype_dict[pop_i]
                right_Gt_i = right_genotype_dict[pop_i]
                sums[:-1, idx] = compute_numer_one_pop_between(
                    left_Gt_i, right_Gt_i, 
                    left_rec_map, right_rec_map, bins,
                    left_mut_map, right_mut_map, u_bar)
            else:
                if not cross_pop:
                    continue
                left_Gt_i = left_genotype_dict[pop_i]
                left_Gt_j = left_genotype_dict[pop_j]
                right_Gt_i = right_genotype_dict[pop_i]
                right_Gt_j = right_genotype_dict[pop_j]
                sums[:-1, idx] = compute_numer_cross_pop_between(
                    left_Gt_i, left_Gt_j, right_Gt_i, right_Gt_j, 
                    left_rec_map, right_rec_map, bins,
                    left_mut_map, right_mut_map, u_bar)
            idx += 1
    sums[-1] = 0
    return sums


def compute_numer_one_pop_within(
    genotypes, 
    rec_map,
    bins, 
    mut_map, 
    u_bar
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
    weights = 1.0 * (genotypes[:, 0] != genotypes[:, 1])
    weights *= (u_bar / mut_map)
    stats = dpluspy.parsing._count_locus_pairs(
        rec_map, 
        bins,
        weights=weights,
        verbose=False
    )
    return stats


def compute_numer_one_pop_between(
    left_genotypes,
    right_genotypes,
    left_rec_map,
    right_rec_map,
    bins, 
    left_mut_map, 
    right_mut_map,
    u_bar
):
    """
    Compute the numerator of the u-adjusted statistic between two genomic
    intervals. This is

    u_bar ** 2 * sum_(i,j) D+(i,j) / (u_i * u_j),

    and the whole estimator is

    u_bar ** 2 / n_pairs * sum_(i,j) D+(i,j) / (u_i * u_j)

    :param genotypes: Array of genotypes for a single diploid
    :param rec_map: Array of recombination map coordinates for genotyped sites
    :param bins: Array of recombination bin edges
    :param mut_map: Array of estimated mutation rates at genotyped sites
    """
    # indicators for one-locus heterozygosity
    left_weights = 1.0 * (left_genotypes[:, 0] != left_genotypes[:, 1])
    left_weights *= (u_bar / left_mut_map)
    right_weights = 1.0 * (right_genotypes[:, 0] != right_genotypes[:, 1])
    right_weights *= (u_bar / right_mut_map)
    stats = dpluspy.parsing._count_locus_pairs_between(
        left_rec_map, 
        right_rec_map,
        bins,
        left_weights=left_weights,
        right_weights=right_weights,
        verbose=False
    )
    return stats


def compute_numer_cross_pop_within(
    genotypes_0, 
    genotypes_1,
    rec_map,
    bins,
    mut_map,
    u_bar
):
    """
    
    """
    weights = dpluspy.parsing._compute_pi_xy(genotypes_0, genotypes_1)
    weights *= (u_bar / mut_map)
    stats = dpluspy.parsing._count_locus_pairs(
        rec_map, 
        bins,
        weights=weights,
        verbose=False
    )
    return stats


def compute_numer_cross_pop_between(
    left_genotypes_0, 
    left_genotypes_1,
    right_genotypes_0,
    right_genotypes_1,
    left_rec_map,
    right_rec_map,
    bins,
    left_mut_map,
    right_mut_map,
    u_bar,
):
    """
    
    """
    left_weights = dpluspy.parsing._compute_pi_xy(
        left_genotypes_0, left_genotypes_1)
    left_weights *= (u_bar / left_mut_map)
    right_weights = dpluspy.parsing._compute_pi_xy(
        right_genotypes_0, right_genotypes_1)
    right_weights *= (u_bar / right_mut_map)
    stats = dpluspy.parsing._count_locus_pairs_between(
        left_rec_map, 
        right_rec_map,
        bins,
        left_weights=left_weights,
        right_weights=right_weights,
        verbose=False
    )
    return stats


def get_genotype_dict(genotypes, sample_ids, sample_labels=None):
    """
    sample_labels[label] = sample_id
    """
    if sample_labels is None:
        sample_mapping = {sample_id: sample_id for sample_id in sample_ids}
    else:
        sample_mapping = sample_labels
    n_samples = genotypes.shape[1]
    assert n_samples == len(sample_ids)
    genotype_dict = dict()
    for sample_label in sample_mapping:
        sample_id = sample_mapping[sample_label]
        idx = sample_ids.index(sample_id)
        genotype_dict[sample_label] = genotypes[:, idx]
    return genotype_dict 


def load_label_file(pop_file):
    """
    Load a population file.
    """
    sample_labels = dict()
    with open(pop_file, 'r') as fin:
        for line in fin:
            label, sample_id = line.split()
            if label in sample_labels:
                raise ValueError("Repeated labels in label file")
            sample_labels[label] = sample_id
    return sample_labels


if __name__ == "__main__": 
    main()
