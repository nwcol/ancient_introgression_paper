"""
Take a site-resolution mutation map (.npy format) and build a windowed analog.
Optionally, apply a genetic mask to the mutation map.
"""

import argparse
import numpy as np
import pandas

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
        help="Pathname of input .npy file")
    parser.add_argument("-scale", "--scale", type=int, required=True,
        help="Window scale, in bp")
    parser.add_argument("-b", "--maskfile", required=True, 
        help="Pathname of optional BED genetic mask file")
    parser.add_argument("-o", "--outfile", required=True, 
        help="Pathname of output .bedgraph file")
    return parser.parse_args()


def build_windowed_map(infile, scale, maskfile, outfile):
    mut_map = np.load(infile)
    L = len(mut_map)
    regions, chrom = utils._read_bed_file(maskfile)
    mask = utils._regions_to_mask(regions, length=L)
    mut_map[mask] = np.nan
    chrom_start = np.arange(0, L, scale)
    chrom_end = np.arange(scale, L + scale, scale)
    windows = np.stack((chrom_start, chrom_end), axis=1)
    windowed_map = np.zeros(len(windows), dtype=np.float64)
    num_sites = np.zeros(len(windows), dtype=np.int64)
    for ii, (start, end) in enumerate(windows):
        # This check exists to avoid raising numpy warnings
        if np.all(np.isnan(mut_map[start:end])):
            windowed_map[ii] = np.nan 
            num_sites[ii] = 0
        else:
            windowed_map[ii] = np.nanmean(mut_map[start:end])
            num_sites[ii] = np.count_nonzero(~np.isnan(mut_map[start:end]))
    data = {
        "chrom": [chrom] * len(windows),
        "chromStart": chrom_start, 
        "chromEnd": chrom_end,
        "num_sites": num_sites, 
        "mut_map": windowed_map
    }
    pandas.DataFrame(data).to_csv(outfile, index=False, na_rep="nan")
    return


def main():
    args = get_args()
    build_windowed_map(args.infile, args.scale, args.maskfile, args.outfile)
    return


if __name__ == "__main__":
    main()
