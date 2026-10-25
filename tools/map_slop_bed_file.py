"""
Extend the intervals in a BED file by some constant map distance, using an
estimate of the recombination map. 
"""

import argparse
import h2py
from h2py.utils import timestamp
import numpy as np
import scipy


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, 
        help='Pathname of input BED file')
    parser.add_argument('-r', '--mapfile', required=True,
        help='Pathname of input recombination map file')
    parser.add_argument('-buffer', '--buffer', required=True, type=float,
        help='Slop distance, in Morgans ~ r')
    parser.add_argument('-o', '--outfile', required=True,
        help='Pathname of output BED file')
    return parser.parse_args()


def slop_bed_file_with_map(infile, mapfile, buffer, outfile):
    """
    Extend BED intervals by `args.buffer` M by loading them, transforming
    them to map coordinates to add/subtract the distance, and transforming them 
    back into physical coordinates.
    """
    intervals, chrom = h2py.masks._read_bed_file(infile)
    pos, coords = h2py.utils._read_rec_map_file(
        mapfile, map_col="cM", pos_col="pos")
    rec_map = scipy.interpolate.interp1d(
        pos, coords, bounds_error=False, fill_value="extrapolate")
    inverse_map = scipy.interpolate.interp1d(
        coords, pos, bounds_error=False, fill_value="extrapolate")
    edge_map = rec_map(intervals)
    edge_map[:, 0] -= buffer
    edge_map[:, 1] += buffer
    slopped = inverse_map(edge_map)
    slopped[:, 0] = np.floor(slopped[:, 0])
    slopped[:, 1] = np.ceil(slopped[:, 1])
    slopped = slopped.astype(np.int64)
    slopped[slopped < 0] = 0
    print(timestamp(), f"chromosome {chrom}")
    print(timestamp(), f"mean slop (lower): {np.mean(intervals[:, 0] - slopped[:, 0])}")
    print(timestamp(), f"mean slop (lower): {np.mean(slopped[:, 1] - intervals[:, 1])}")
    # Remove any redundant regions
    slopped_mask = h2py.GeneticMask.from_intervals(slopped, chrom=chrom)
    slopped_mask.to_bed_file(outfile)
    return


def main():
    args = get_args()
    slop_bed_file_with_map(
        args.infile, args.mapfile, args.buffer, args.outfile)
    return


if __name__ == '__main__':
    main()
