## utilities, mostly for reading and writing common types of file

from datetime import datetime
import gzip
import numpy as np
import re
import warnings


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


def intersect_bed_regions(*bed_regions):
    """

    """
    length = max([reg[-1, 1] for reg in bed_regions])
    masks = [regions_to_mask(reg, length=length) for reg in bed_regions]
    sums = np.sum(masks, axis=1)
    overlap_mask = sums < 0
    regions = mask_to_regions(overlap_mask)
    return regions


def collapse_regions(elements):
    """
    Collapse any overlapping elements in an array together.
    """
    return mask_to_regions(regions_to_mask(elements))


def read_bedgraph(file, sep=','):
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
    # check for proper header format
    assert header_line[0] in ['chrom', '#chrom']
    assert header_line[1] in ['chromStart', 'start']
    assert header_line[2] in ['chromEnd', 'end']
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


def write_bedgraph(file, chrom_num, regions, data, sep=','):
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
