"""
Extend the intervals in a BED file by some constant map distance, using an
estimate of the recombination map. 
"""

import argparse
import numpy as np

from dpluspy import utils, parsing


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', required=True, 
        help='Pathname of input BED file')
    parser.add_argument('-r', '--map_file', required=True,
        help='Pathname of input recombination map file')
    parser.add_argument('-buffer', '--buffer', required=True, type=float,
        help='Slop distance, in Morgans ~ r')
    parser.add_argument('-o', '--out_file', required=True,
        help='Pathname of output BED file')
    return parser.parse_args()


def main():
    """
    Extend BED intervals by `args.buffer` M by loading them, transforming
    them to map coordinates to add/subtract the distance, and transforming them 
    back into physical coordinates.
    """
    args = get_args()
    regions, chrom_num = utils._read_bed_file(args.in_file)
    rec_map = parsing._load_recombination_map(
        args.map_file, map_col="cM", pos_col="pos")
    inverse_map = parsing._load_recombination_map(
        args.map_file, inverse=True, map_col="cM", pos_col="pos")
    buffer = args.buffer
    edge_map = rec_map(regions)
    edge_map[:, 0] -= buffer
    edge_map[:, 1] += buffer 
    slopped = inverse_map(edge_map)
    slopped[:, 0] = np.floor(slopped[:, 0])
    slopped[:, 1] = np.ceil(slopped[:, 1])  
    slopped = slopped.astype(np.int64)
    slopped[slopped < 0] = 0
    # Remove any redundant regions
    slopped = utils._mask_to_regions(utils._regions_to_mask(slopped))
    utils._write_bed_file(args.out_file, slopped, chrom_num)

    return


if __name__ == '__main__':
    main()
