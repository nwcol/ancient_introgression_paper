## utilities, mostly for reading and writing common types of file

from datetime import datetime
import gzip
import numpy as np
import re
import pickle
import warnings




#### NEW STUFF



## Generating names of statistics


def generate_pairs(pop_ids):

    pairs = []
    for i, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[i:]:
            pairs.append((pop_i, pop_j))

    return pairs


def h_names(pop_ids):

    hs = []
    for i, pop_id0 in enumerate(pop_ids):
        for pop_id1 in pop_ids[i:]:
            hs.append(f"H_{pop_id0}_{pop_id1}")

    return hs


def d_plus_names(pop_ids):

    ds = []
    for i, pop_id0 in enumerate(pop_ids):
        for pop_id1 in pop_ids[i:]:
            ds.append(f"D+_{pop_id0}_{pop_id1}")

    return ds


def stat_names(pop_ids):

    ds = d_plus_names(pop_ids)
    hs = h_names(pop_ids)

    return (ds, hs)


def h_idxs(num_pops):

    hs = []
    for i in range(num_pops):
        for j in range(i, num_pops):
            hs.append(f"H_{i}_{j}")

    return hs


def d_plus_idxs(num_pops):

    ds = []
    for i in range(num_pops):
        for j in range(i, num_pops):
            ds.append(f"D+_{i}_{j}")

    return ds


def get_latex_names(pop_ids, statistic="D^+"):
    """
    From a list of population names, get a list of strings of the form 
    '$statistic_{pop0,pop1}$' for each pair of populations.

    :param pop_ids: List of population names.
    :type pop_ids: list
    :param statistic: Name of the statistic (default 'D^+')
    :type statistic: str

    :returns: A list of string statistic names in a LaTeX-friendly format.
    :rtype: list
    """
    names = []
    for i, pop0 in enumerate(pop_ids):
        for pop1 in pop_ids[i:]:
            if pop0 == pop1:
                names.append(f"${statistic}_{{{pop0}}}$")
            else:
                names.append(f"${statistic}_{{{pop0},{pop1}}}$")

    return names


## Subsetting empirical statistics


def subset_oldstyle_stats(statistics, to_pops, min_r=None, max_r=None):
    # for old format. stats is a dictionary
    means = statistics["means"]
    varcovs = statistics["covs"]
    pop_ids = statistics["pops"]
    if min_r is not None or max_r is not None:
        raise ValueError("not implemented")
    new_means = subset_means(means, pop_ids, to_pops)
    new_varcovs = subset_varcovs(varcovs, pop_ids, to_pops)
    bins = statistics["bins"]

    return bins, new_means, new_varcovs


def load_statistics(filename, to_pops=None):
    
    stats = pickle.load(open(filename, "rb"))
    bins = stats["bins"]
    if to_pops is not None:
        all_ids = stats["pop_ids"]
        means = subset_means(stats["means"], all_ids, to_pops)
        varcovs = subset_varcovs(stats["varcovs"], all_ids, to_pops)
        pop_ids = to_pops
    else:
        pop_ids = stats["pop_ids"]
        means = stats["means"]
        varcovs = stats["varcovs"]

    return pop_ids, bins, means, varcovs


def subset_means(means, pop_ids, to_pops):
    """
    Marginalize statistics for `pop_ids` to `pops`. 

    :param means: List of 1d arrays to subset.  
    :type means: list
    :param pop_ids: List of populations represented in `means`.
    :type pop_ids: list
    :param to_pops: List of populations to subset to. One and two-population
        statistics from this list will be returned.
    :typr to_pops: list

    :returns: A list of 1d arrays subset to `to_pops`.
    """
    for pop in to_pops:
        if pop not in pop_ids:
            raise ValueError(f"{pop} not in `pop_ids`")
    stats = d_plus_idxs(len(pop_ids))
    to_pop_idx = [pop_ids.index(pop) for pop in to_pops]
    to_stats = []
    for i, idx0 in enumerate(to_pop_idx):
        for idx1 in to_pop_idx[i:]:
            _idx0, _idx1 = sorted([idx0, idx1])
            to_stats.append(f"D+_{_idx0}_{_idx1}")
    to_idx = np.array([stats.index(to_stat) for to_stat in to_stats])
    new_means = [means[i][to_idx] for i in range(len(means))]

    return new_means


def subset_varcovs(varcovs, pop_ids, to_pops):
    """
    Marginalize covariance matrices from `pop_ids` to `pops`.

    :returns: A list of 2d covariance matrices subset to `to_pops`.
    """
    for pop in to_pops:
        if pop not in pop_ids:
            raise ValueError(f"{pop} not in `pop_ids`")
    stats = d_plus_idxs(len(pop_ids))
    to_pop_idx = [pop_ids.index(pop) for pop in to_pops]
    to_stats = []
    for i, idx0 in enumerate(to_pop_idx):
        for idx1 in to_pop_idx[i:]:
            _idx0, _idx1 = sorted([idx0, idx1])
            to_stats.append(f"D+_{_idx0}_{_idx1}")
    to_idx = np.array([stats.index(to_stat) for to_stat in to_stats])
    mesh = np.ix_(to_idx, to_idx)
    new_varcovs = [varcovs[i][mesh] for i in range(len(varcovs))]

    return new_varcovs












### functions moved here from parsing module


def subset_statistics(
    data, 
    graph=None, 
    to_pops=None, 
    min_dist=None, 
    max_dist=None
):
    """
    Subset a dictionary of statistics by `pop`. If a graph is provided, subsets
    to the set of names which occur in both the data set and the graph.
    """
    ret = copy.deepcopy(data)
    
    if graph is not None:
        if to_pops is not None:
            warnings.warn('argument `to_pops` overriden by `graph`')
        if isinstance(graph, str):
            graph = demes.load(graph)
        graph_demes = [d.name for d in graph.demes]
        pops = data['pops']
        to_pops = [d for d in pops if d in graph_demes]

    if to_pops is not None:
        pops = data['pops']

        # sort `to_pops` so it's in the same order as `pops`
        indices = [pops.index(p) for p in to_pops]
        to_pops = [pops[i] for i in sorted(indices)]

        two_pop = False if data['means'].shape[1] == len(pops) else True
        labels = enumerate_labels(pops=pops, two_pop=two_pop)
        to_labels = enumerate_labels(pops=to_pops, two_pop=two_pop)
        keep = np.array([labels.index(label) for label in to_labels])

        for key in ['sums', 'means']:
            if key in data:
                ret[key] = data[key][:, keep]
        
        if 'covs' in data and data['covs'] is not None:
            covs = data['covs']
            ret['covs'] = np.stack([cov[np.ix_(keep, keep)] for cov in covs])
        
        ret['pops'] = to_pops

    if min_dist is not None or max_dist is not None:
        bins = data['bins']
        if min_dist is None:
            min_dist = 0
        if max_dist is None:
            max_dist = np.inf
        min_bin = np.searchsorted(bins, min_dist)
        max_bin = np.searchsorted(bins, max_dist)
        
        data['bins'] = bins[min_bin:max_bin + 1]

        for key in ['sums', 'means', 'denom']:
            if key in data:
                ret[key] = ret[key][min_bin:max_bin + 1]
        
        if 'covs' in data:
            ret['covs'] = ret['covs'][min_bin:max_bin + 1]

    return ret


















## .bed files and genetic masks


def read_bedfile(bed_file, return_chrom=False):
    """
    Load the regions in a .bed file as a (num_regions, 2) shape numpy array.
    """
    open_func = gzip.open if bed_file.endswith('.gz') else open

    with open_func(bed_file, 'rb') as fin:
        split_line = fin.readline().decode().split()
        if split_line[1].isnumeric():
            skiprows = 0
        else:
            skiprows = 1

        if return_chrom:
            if skiprows == 0:
                chrom_num = split_line[0]
            else:
                chrom_num = fin.readline().decode().split()[0]

    regions = np.loadtxt(bed_file, usecols=(1, 2), dtype=int, skiprows=skiprows)

    if regions.ndim == 1:
        regions = regions[np.newaxis, :]

    if return_chrom:
        ret = (regions, chrom_num)
    else:
        ret = regions

    return ret


def read_bedfile_positions(bed_file, return_chrom=False):
    """
    Get the vector of 0-indexed positions within the regions specified in 
    `bed_file`.
    """
    regions = read_bedfile(bed_file, return_chrom=False)
    mask = regions_to_mask(regions)
    positions = np.nonzero(~mask)[0]
    
    return positions


def write_bedfile(file, regions, chrom_num=None, chrom_nums=None, header=False):
    """
    
    """
    open_func = gzip.open if file.endswith('.gz') else open

    if chrom_num is not None:
        chr_col = [chrom_num] * len(regions)

    elif chrom_nums is not None:
        assert len(chrom_nums) == len(regions)
        chr_col = chrom_nums

    else:
        raise ValueError('please provide chromosome numbers')

    with open_func(file, "wb") as file:
        if header:
            header = b'#chrom\tchromStart\tchromEnd\n'
            file.write(header)
            
        for i, (start, stop) in enumerate(regions):
            line = f'{chr_col[i]}\t{start}\t{stop}\n'.encode()
            file.write(line)
    return 


def regions_to_mask(regions, length=None):
    """
    Return a boolean mask array that equals 0 within `regions` and 1 
    elsewhere.
    """
    if length is None:
        length = regions[-1, 1]
    mask = np.ones(length, dtype=bool)
    for (start, end) in regions:
        if start >= length:
            continue
        elif end > length:
            end = length
        mask[start:end] = 0
    return mask


def mask_to_regions(mask):
    """
    Return an array representing the regions that are not masked in a boolean
    array (0s).
    """
    jumps = np.diff(np.concatenate(([1], mask, [1])))
    starts = np.where(jumps == -1)[0]
    ends = np.where(jumps == 1)[0]
    regions = np.stack([starts, ends], axis=1)
    return regions


def intersect_regions(*regionss):
    """

    """
    length = max([reg[-1, 1] for reg in regionss])
    masks = [regions_to_mask(reg, length=length) for reg in regionss]
    sums = np.sum(masks, axis=1)
    overlap_mask = sums < 0
    regions = mask_to_regions(overlap_mask)
    return regions


def collapse_regions(elements):
    """
    Collapse any overlapping elements in an array together.
    """
    return mask_to_regions(regions_to_mask(elements))


def read_bedgraph(file, sep='\t'):
    """
    From a bedgraph-format file, read and return chromosome number(s), an 
    array of genomic regions and a dictionary of data columns. 

    If the file has one unique chromosome number, returns it as a string of
    the form `chr00`; if there are several, returns an array of string
    chromosome numbers of this form for each row.
    Possible file extensions include but are not limited to .bedgraph, .csv,
    and .tsv, with column seperator determined by the `sep` argument.
    """
    open_func = gzip.open if file.endswith('.gz') else open
    with open_func(file, 'rb') as fin:
        header_line = fin.readline().decode().strip().split(sep)
    fields = header_line[3:]
    # handle the return of the chromosome number(s)
    chrom_nums = np.loadtxt(
        file, usecols=0, dtype=str, skiprows=1, delimiter=sep
    )
    if len(set(chrom_nums)) == 1:
        ret_chrom = chrom_nums[0]
    else:
        # return the whole vector if there are >1 unique chromosome
        ret_chrom = chrom_nums
    windows = np.loadtxt(
        file, usecols=(1, 2), dtype=int, skiprows=1, delimiter=sep
    )
    cols_to_load = tuple(range(3, len(header_line)))
    arr = np.loadtxt(
        file,
        usecols=cols_to_load,
        dtype=float,
        skiprows=1,
        unpack=True,
        delimiter=sep
    )
    dataT = [arr] if arr.ndim == 1 else [col for col in arr]
    data = dict(zip(fields, dataT))
    return windows, data


def write_bedgraph(file, chrom_num, regions, data, sep='\t'):
    """
    Write a .bedgraph-format file from an array of regions/windows and a 
    dictionary of data columns.
    """
    for field in data:
        if len(data[field]) != len(regions):
            raise ValueError(f'data field {data} mismatches region length!')
    open_func = gzip.open if file.endswith('.gz') else open
    fields = list(data.keys())
    header = sep.join(['#chrom', 'chromStart', 'chromEnd'] + fields) + '\n'
    with open_func(file, 'wb') as file:
        file.write(header.encode())
        for i, (start, end) in enumerate(regions):
            ldata = [str(data[field][i]) for field in fields]
            line = sep.join([chrom_num, str(start), str(end)] + ldata) + '\n'
            file.write(line.encode())

    return


## dealing with recombination maps


def get_uniform_rec_map(r, sites):
    """
    Obtain a recombination map for `sites` assuming a constant recombination
    rate `r`. Returns a map in units of cM.
    """
    cM_per_bp = map_function(r)
    rec_map = sites * cM_per_bp

    return rec_map


def get_rec_map(rec_map_file, sites, map_col="Map(cM)"):
    
    if (
        rec_map_file.endswith(".txt") 
        or rec_map_file.endswith(".txt.gz")
    ):
        coords, vals = read_hapmap_rec_map(rec_map_file, map_col=map_col)

    elif (
        rec_map_file.endswith(".bedgraph") 
        or rec_map_file.endswith(".bedgraph.gz")
    ):
        coords, vals = read_bedgraph_rec_map(rec_map_file)

    else:
        raise ValueError("unrecognized recombination map file type")

    assert np.all(np.diff(coords) > 0)
    assert np.all(np.diff(vals) >= 0)
    site_map = np.interp(sites, coords, vals)

    return site_map


def read_recombination_map(rec_map_file, positions, map_col='Map(cM)'):
    """
    Read a recombination map in `hapmap` map format and interpolate map values
    for a vector of positions.
    """
    open_func = gzip.open if rec_map_file.endswith('.gz') else open
    with open_func(rec_map_file, 'rb') as fin:
        header_line = fin.readline().decode().split()
    pos_idx = header_line.index('Position(bp)')
    map_idx = header_line.index(map_col)
    map_coords, map_vals = np.loadtxt(
        rec_map_file, skiprows=1, usecols=(pos_idx, map_idx), unpack=True
    )

    assert np.all(np.diff(map_coords) > 0)
    assert np.all(np.diff(map_vals) >= 0)
    
    pos_map = np.interp(positions, map_coords, map_vals)
    return pos_map


def read_hapmap_rec_map(map_file, map_col="Map(cM)"):
    ## get coords and map values
    open_func = gzip.open if map_file.endswith('.gz') else open

    with open_func(map_file, 'rb') as fin:
        header_line = fin.readline().decode().split()

    coord_idx = header_line.index('Position(bp)')
    map_idx = header_line.index(map_col)

    coords = np.loadtxt(
        map_file, skiprows=1, usecols=(coord_idx), dtype=np.int64
    )
    vals = np.loadtxt(map_file, skiprows=1, usecols=(map_idx), dtype=np.float64)

    return coords, vals


def read_bedgraph_rec_map(map_file):
    # assumes map coordinate is in the LAST column, e.g.
    # chrom start end ... map_coord
    open_func = gzip.open if map_file.endswith('.gz') else open

    coords = []
    mapvals = []

    with open_func(map_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith("#"):
                continue
        
            split_line = line.split()
            coords.append(int(split_line[2]))
            mapvals.append(float(split_line[-1]))

    coords = np.array(coords, dtype=np.int64)
    mapvals = np.array(mapvals, dtype=np.float64)

    return coords, mapvals


def read_mutation_map(mut_map_file, positions):

    """
    """
    if (
        mut_map_file.endswith('.bedgraph') 
        or mut_map_file.endswith('.bedgraph.gz')
    ):
        regions, data = read_bedgraph(mut_map_file)
        # interpolate.
        idxs = np.searchsorted(regions[:, 1], positions)
        reg_mut_map = data['mut_rate']
        mut_map = reg_mut_map[idxs]
        
    elif mut_map_file.endswith('.npy'):
        tot_mut_map = np.load(mut_map_file)
        mut_map = tot_mut_map[positions]
        assert not np.any(np.isnan(mut_map))

    else:
        raise ValueError('unrecognized mutation map format')

    return mut_map


## reading .vcf files (under construction)


def extract_gqs(samples, gq_index):
    """
    
    """
    gq_strs = [s.split(':')[gq_index] for s in samples]
    gqs = np.array([np.nan if gq == '.' else float(gq) for gq in gq_strs])
    return gqs


def extract_genotypes(samples, missing_to_ref=True):
    """
    not finished
    """
    if missing_to_ref:
        gts = [re.split('/|\|', s.split(':')[0]) for s in samples]
        if '.' in str(gts):
            gts = [['0', '0'] if '.' in x else x for x in gts]
            print('replaced missing alleles')
        genotypes = np.array(gts, dtype=np.int64)
    else:
        gts = [re.split('/|\|', s.split(':')[0]) for s in samples]
        genotypes = np.array(gts, dtype=np.int64)
    return genotypes


def read_genotypes(
    vcf_file, 
    bed_file=None, 
    min_reg_len=None,
    region=None,
    ancestral_seq=None, 
    read_multiallelic=False,
    missing_to_ref=True,
):
    """
    Read a genotype matrix from a .vcf file. Matrix has shape 
    (nsamples, nsites, 2). Ignores sites that are not biallelic. 
    """
    if bed_file is not None:
        _, regions = read_bedfile(bed_file)
        mask = regions_to_mask(regions)
        len_mask = len(mask)
    
    open_func = gzip.open if vcf_file.endswith('.gz') else open

    first_row = True
    outside_mask = 0
    outside_reg = 0
    multiallelic = 0

    chrom_nums = []
    positions = []
    genotypes = []

    with open_func(vcf_file, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                if line.startswith('#CHROM'):
                    sample_ids = line.split()[9:]
                    num_samples = len(sample_ids)
            else:
                split_line = line.split()
                chrom, pos, _, ref, alt = split_line[:5]

                if first_row:
                    fmt = split_line[8]
                    if 'GQ' in fmt:
                        gq_index = fmt.split(':').index('GQ')
                    else:
                        gq_index = None
                        
                    first_row = False

                samples = split_line[9:]
                position = int(pos) - 1

                if region is not None:
                    if position < region[0] or position >= region[-1]:
                        outside_reg += 1
                        continue

                if bed_file is not None:
                    if position >= len_mask or mask[position] == 1:
                        outside_mask += 1
                        continue

                if not read_multiallelic:
                    if len(alt.split(',')) > 1 or len(ref) > 1:
                        multiallelic += 1
                        continue

                line_genotypes = extract_genotypes(
                    samples, missing_to_ref=missing_to_ref
                )
                chrom_nums.append(chrom)
                positions.append(position)
                genotypes.append(line_genotypes)

    positions = np.array(positions, dtype=np.int64)
    genotypes = np.stack(genotypes, axis=1, dtype=np.int64)
    unique_chrom_nums = list(set(chrom_nums))
    if len(unique_chrom_nums) > 1:
        warnings.warn('more than one unique chromosome in .vcf')
    chrom_num = unique_chrom_nums[0]

    return chrom_num, sample_ids, positions, genotypes


## math


def n_choose_2(n):
    """
    
    """
    return n * (n - 1) // 2


## recombination map math


def map_function(r):
    """
    Haldane's map function; transforms distance in r to cM.
    """
    return -50 * np.log(1 - 2 * r)


def inverse_map_function(d):
    """
    The inverse of Haldane's map function. Transforms distance in cM to r.
    """
    return (1 - np.exp(-d / 50)) / 2


## printouts


def get_time():
    """
    Return a string giving the time and date with yy-mm-dd format.
    """
    return "[" + datetime.strftime(datetime.now(), "%y-%m-%d %H:%M:%S") + "]"


## msprime simulation helper functions


def increment1(x):
    """
    Increment 1 to every .vcf site to make it 1-indexed rather than 0-indexed.
    """
    return [_ + 1 for _ in x]

