"""
Extend the intervals in a BED file by some constant map distance, using an
estimate of the recombination map. 
"""

import argparse
import dpluspy
import numpy as np


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
    args = get_args()
    regions, chrom_num = dpluspy.utils._read_bed_file(infile)
    rec_map = dpluspy.parsing._load_recombination_map(
        mapfile, map_col="cM", pos_col="pos")
    inverse_map = dpluspy.parsing._load_recombination_map(
        mapfile, inverse=True, map_col="cM", pos_col="pos")
    buffer = buffer
    edge_map = rec_map(regions)
    edge_map[:, 0] -= buffer
    edge_map[:, 1] += buffer 
    slopped = inverse_map(edge_map)
    slopped[:, 0] = np.floor(slopped[:, 0])
    slopped[:, 1] = np.ceil(slopped[:, 1])  
    slopped = slopped.astype(np.int64)
    slopped[slopped < 0] = 0
    # Remove any redundant regions
    slopped = dpluspy.utils._mask_to_regions(
        dpluspy.utils._regions_to_mask(slopped))
    dpluspy.utils._write_bed_file(outfile, slopped, chrom_num)
    return


def main():
    args = get_args()
    slop_bed_file_with_map(
        args.infile, args.mapfile, args.buffer, args.outfile)
    return


if __name__ == '__main__':
    main()
