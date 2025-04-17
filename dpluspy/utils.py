"""
Utilities for reading/writing files and manipulating statistics, arrays
"""

import copy
from datetime import datetime
import gzip
import numpy as np
import re
import pickle
import warnings


## Generating the names of statistics


def generate_pairs(pop_ids):
    """
    Generate a list of 2-tuples holding the (n choose 2) unique pairs that may
    be drawn from `pop_ids`. 

    :dtype pop_ids: list
    :rtype: list of tuples
    """
    pairs = []
    for i, pop_i in enumerate(pop_ids):
        for pop_j in pop_ids[i:]:
            pairs.append((pop_i, pop_j))

    return pairs


def H_names(pop_ids):
    """
    Generate a list of names of the unique one and two-population H statistics 
    corresponding to a list of population IDs. 

    :dtype pop_ids: list of strings
    :rtype: list of strings
    """
    names = []
    for i, pop_id0 in enumerate(pop_ids):
        for pop_id1 in pop_ids[i:]:
            names.append(f"H_{pop_id0}_{pop_id1}")

    return names


def Dplus_names(pop_ids):
    """
    Get a list of the unique one and two-population D+ statistics corresponding
    to a list of population IDs.

    :dtype pop_ids: list of strings
    :rtype: list of strings
    """
    names = []
    for i, pop_id0 in enumerate(pop_ids):
        for pop_id1 in pop_ids[i:]:
            names.append(f"D+_{pop_id0}_{pop_id1}")

    return names


def stat_names(pop_ids):
    """
    Get the names of all the D+ and H statistics for populations `pop_ids`.
    Statistic names have the form 'D+_{pop_i}_{pop_j}' and 'H_{pop_i}_{pop_j}'.

    :param pop_ids: List of population names.
    :type pop_ids: list of str

    :returns: Lists of names for D+ and H statistics.
    :rtype: tuple of lists of strings
    """
    Dplus_names = Dplus_names(pop_ids)
    H_names = H_names(pop_ids)

    return (Dplus_names, H_names)


def get_latex_names(pop_ids, statistic="D^+"):
    """
    From a list of population names, get a list of strings of the form 
    '${statistic}_{pop0,pop1}$' for each pair of populations.

    :param pop_ids: List of population names.
    :type pop_ids: list
    :param statistic: Name of the statistic (default 'D^+')
    :type statistic: str

    :returns: A list of string statistic names in a LaTeX-friendly format.
    :rtype: list of strings
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


def subset_statistics(
    statistics, 
    to_pops=None, 
    min_r=None, 
    max_r=None,
    return_dict=False
):
    """
    Subset a dictionary holding statistics by populations or bins. 

    :param statistics: A dictionary with fields 'means', 'varcovs', 'pop_ids',
        and 'bins'.
    :type statistics: dict
    :param to_pops: List of population IDs to subset to (default None).
    :param min_r: Minimum lower bin edge, inclusive (default None).
    :param max_r: Maximum upper bin edge, inclusive (default None).
    :param return_dict: If True, return a dictionary with the same fields as
        required for the input- otherwise return bins, means and varcovs in a 
        tuple (default False).

    :returns: Dictionary of subsetted statistics.
    :rtype: dict
    """
    means = statistics['means']
    varcovs = statistics['varcovs']
    pop_ids = statistics['pop_ids']
    if to_pops is None:
        to_pops = pop_ids
    for pop_id in to_pops:
        if pop_id not in pop_ids:
            raise ValueError(f'"{pop_id}" is not represented in the data')
    if min_r is not None or max_r is not None:
        if min_r is not None:
            min_idx = np.where(bins >= min_r)[0][0]
        else:
            min_idx = 0
        if max_r is not None:
            max_idx = np.where(bins <= max_r)[0][-1]
        else:
            max_idx = len(bins) - 1
        means = means[min_idx:max_idx] + [means[-1]]
        varcovs = varcovs[min_idx:max_idx] + [varcovs[-1]]
    else:
        bins = statistics['bins']
    new_means = subset_means(means, pop_ids, to_pops)
    new_varcovs = subset_varcovs(varcovs, pop_ids, to_pops)

    if return_dict:
        subset_stats = {
            'pop_ids': pop_ids,
            'bins': bins,
            'means': new_means,
            'varcovs': new_varcovs
        }
        return subset_stats
    else:
        return bins, new_means, new_varcovs


def load_statistics(filename, to_pops=None):
    # deprecated
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
    Subset a list of binned statistics representing `pop_ids` to `to_pops`. 

    :param means: List of 1d arrays to subset.  
    :type means: list of np.ndarray
    :param pop_ids: List of populations represented in `means`.
    :type pop_ids: list of str
    :param to_pops: List of populations to subset to. One and two-population
        statistics from this list will be returned.
    :typr to_pops: list of str

    :returns: A list of 1d arrays subset to `to_pops`.
    :rtype: list of np.ndarray
    """
    for pop in to_pops:
        if pop not in pop_ids:
            raise ValueError(f'"{pop}" not in `pop_ids`')
    stats = Dplus_names(len(pop_ids))
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
    Marginalize a list of bin-wise covariance matrices from `pop_ids` to `pops`.

    :returns: A list of 2d covariance matrices subset to `to_pops`.
    :rtype: list of np.ndarray
    """
    for pop in to_pops:
        if pop not in pop_ids:
            raise ValueError(f'"{pop}" not in `pop_ids`')
    stats = Dplus_names(len(pop_ids))
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


## BED files and genetic masks


def read_bed_file(filename):
    """
    Load regions from a BED file as an array. Expects the structure 
        CHROM\tSTART\tEND...\n
    on each line, and skips any comment/header lines that begin with '#'. 
    Raises an error if the BED file has more than one unique CHROM entry.

    :param filename: Pathname of the BED file to load. 
    :type filename: str

    :returns: ndarray of BED file regions, BED chromosome ID
    :rtype: np.ndarray, str
    """
    if filename.endswith('.gz'):
        openfunc = gzip.open 
    else:
        openfunc = open
    chroms = []
    starts = []
    ends = []
    with openfunc(filename, "rb") as fin:
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                continue
            split_line = line.split()
            chroms.append(split_line[0])
            starts.append(split_line[1])
            ends.append(split_line[2])
    chrom_set = set(chroms)
    # check that there is one unique CHROM
    if len(chrom_set) > 1:
        raise ValueError('BED files must describe one chromosome only')
    # check to make sure one or more lines were read
    elif len(chrom_set) == 0:
        raise ValueError('BED file has no valid contents')
    chrom = list(chrom_set)[0]
    regions = np.array([[start, end] for start, end in zip(starts, ends)])

    return regions, chrom


def read_bed_file_positions(bed_file):
    """
    Read a BED file and return a vector of the positions recorded in its
    intervals (0-indexed).
    """
    regions = read_bed_file(bed_file)[0]
    mask = regions_to_mask(regions)
    positions = np.nonzero(~mask)[0]
    
    return positions


def write_bed_file(filename, regions, chrom):
    """
    Write a BED file. Does not write a header.

    :param filename: Pathname of output file. Should end in .bed or .bed.gz. 
    :type filename: str
    :param regions: Array of BED regions to save.
    :type regions: np.ndarray
    :param chrom: Chromosome number to use in the CHROM column. All regions are 
        assigned to the same chromosome. 
    :type chrom: str

    :returns: None
    """
    if filename.endswith('.gz'):
        openfunc = gzip.open 
    else:
        openfunc = open
    with openfunc(filename, 'wb') as fout:
        fout.write('#CHROM\tSTART\tEND\n'.encode())
        for start, end in regions:
            fout.write(f'{chrom}\t{start}\t{end}\n'.encode())

    return 


def regions_to_mask(regions, length=None):
    """
    Return a boolean mask array that equals False within intervals in `regions` 
    and True elsewhere.

    :param regions: Array of intervals.
    :param length: Optional maximum mask length (default None).

    :returns: Boolean mask array.
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
    Return an array of intervals that equal False in a boolean mask array
    (0-indexed).
    """
    jumps = np.diff(np.concatenate(([1], mask, [1])))
    starts = np.where(jumps == -1)[0]
    ends = np.where(jumps == 1)[0]
    regions = np.stack([starts, ends], axis=1)

    return regions


def intersect_regions(regions_arrs):
    """
    Build an array of intervals where every input regions array has coverage.

    :param regions_arrs: List of BED region arrays.
    :type regions_arrs: list of np.ndarray

    :returns: Region array representing intersection of sites in inputs.
    :rtype: np.ndarray
    """
    length = max([reg[-1, 1] for reg in regions_arrs])
    masks = [regions_to_mask(region, length=length) for region in regions_arrs]
    sums = np.sum(masks, axis=1)
    overlap_mask = sums < 0
    regions = mask_to_regions(overlap_mask)

    return regions


def collapse_regions(regions):
    """
    Collapse any overlapping intervals in an array together.
    """
    return mask_to_regions(regions_to_mask(regions))


# BEDGRAPH files and recombination maps


def read_bedgraph_file(filename, sep=None, override_cols=None):
    """
    From a bedgraph-format file, read and return an array of genomic intervals,
    a dictionary of data and the associated chromosome number. There must be
    a header with format corresponding to
        Chrom\tchromStart\tchromEnd\tdata_col1\t...\n
    in the first line. Other commented or header lines beginning with '#' will
    be ignored.

    :param filename: Pathname of the file to load.
    :param sep: File seperator to expect (default None uses \t).
    :param override_cols: If given, overrides the data field names in the file 
        header (default None).

    :returns: Array of intervals, dictionary of data arrays, and chromosome ID
    """
    if sep is None:
        sep = '\t'
    if filename.endswith('.gz'):
        openfunc = gzip.open 
    else:
        openfunc = open
    chroms = []
    starts = []
    ends = []
    with openfunc(filename, "rb") as fin:
        header_line = fin.readline().decode()
        if header_line[0] != '#':
            raise ValueError('Input file lacks a header line')
        split_header = header_line.strip().split(sep)
        if override_cols is not None:
            if len(override_cols) != len(split_header) - 3:
                raise ValueError('Invalid `override_cols`')
        raw_data = {i: [] for i in range(3, len(split_header))}
        for lineb in fin:
            line = lineb.decode()
            if line.startswith('#'):
                continue
            split_line = line.strip().split(sep)
            for idx in raw_data:
                raw_data[idx].append(split_line[idx])
            chroms.append(split_line[0])
    chrom_set = set(chroms)
    # check that there is one unique CHROM
    if len(chrom_set) > 1:
        raise ValueError('BED files must describe one chromosome only')
    # check to make sure one or more lines were read
    elif len(chrom_set) == 0:
        raise ValueError('BED file has no valid contents')
    chrom = list(chrom_set)[0]
    data = {}
    for idx in raw_data:
        if override_cols is None:
            field = split_header[idx]
        else:
            field = override_cols[idx - 3]
        if '.' in raw_data[idx][0]:
            arr = np.array(raw_data[idx], dtype=np.float64)
        else:
            arr = np.array(raw_data[idx], dtype=np.int64)
        data[field] = arr
    regions = np.array(
        [[start, end] for start, end in zip(starts, ends)], dtype=np.int64
    )
    return regions, data, chrom


def write_bedgraph_file(filename, regions, data, chrom_num, sep=None):
    """
    Write a .bedgraph-format file from an array of regions/windows and a 
    dictionary of data columns.
    """
    for field in data:
        if len(data[field]) != len(regions):
            raise ValueError(f'data field {data} mismatches region length!')
    if sep is None:
        '\t'
    if filename.endswith('.gz'):
        openfunc = gzip.open 
    else:
        openfunc = open
    constants = ['#chrom', 'chromStart', 'chromEnd']
    fields = list(data.keys())
    header = sep.join(constants + fields) + '\n'
    with openfunc(file, 'wb') as file:
        file.write(header.encode())
        for i, (start, end) in enumerate(regions):
            interval = [chrom_num, str(start), str(end)]
            line_data = [str(data[field][i]) for field in fields]
            line = sep.join(interval + line_data) + '\n'
            file.write(line.encode())

    return


def read_bedgraph_map(filename, map_col=None, sep=None):
    """
    Read a map from a BEDGRAPH file, returning an array of physical and of map
    coordinates. If no `map_col` is given, accesses the rightmost column.
    """
    intervals, data = read_bedgraph_file(filename, sep=sep)
    coords = intervals[0, :]
    map_coords = data[map_col]

    return coords, map_coords


def read_hapmap_map(filename, map_col=None):
    """
    Read a recombination map in the Hapmap format, returning arrays of physical
    and map coordinates. The first line must be a header.
    """
    if map_col is None:
        map_col = 'Map(cM)'
    if filename.endswith('.gz'):
        openfunc = gzip.open 
    else:
        openfunc = open 
    coords = []
    map_coords = []
    with openfunc(filename, "rb") as fin:
        header_line = fin.readline().decode().split()
        coord_idx = header_line.index('Position(bp)')
        map_idx = header_line.index(map_col)
        for line in fin:
            split_line = line.decode().strip().split()
            coords.append(split_line[coord_idx])
            map_coords.append(split_line[map_idx])
    coords = np.array(coords, dtype=np.int64)
    map_coords = np.array(map_coords, dtype=np.float64)

    return coords, map_coords


## Recombination map math


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


## Math


def n_choose_2(n):
    """
    Return (n choose 2).
    """
    return n * (n - 1) // 2


## Printouts


def get_time():
    """
    Return a string giving the time and date with yy-mm-dd format.
    """
    return '[' + datetime.strftime(datetime.now(), '%y-%m-%d %H:%M:%S') + ']'


## Simulation helper functions


def increment1(x):
    """
    Increment 1 to every value of x. Used to transform simulated VCF sites to
    1-indexing.
    """
    return [_ + 1 for _ in x]

