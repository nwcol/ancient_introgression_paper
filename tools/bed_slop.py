"""
Extend the regions in a BED file by some distance and save the result.
"""

import argparse
import numpy as np

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-bed', '--bed_file', required=True, 
        help='Pathname of input BED file')
    parser.add_argument('-dist', '--distance', required=True, type=int,
        help='Slop distance, in bp')
    parser.add_argument('-out', '--out_file', required=True,
        help='Pathname of output BED file')
    return parser.parse_args()


def main():
    """
    Load a BED file and extend its intervals by `args.distance`.
    """
    args = get_args()
    regions, chrom_num = utils._read_bed_file(args.bed_file)
    dist = args.distance
    regions[:, 0] -= dist
    regions[:, 1] += dist
    regions[regions < 0] = 0
    resolved = utils._mask_to_regions(utils._regions_to_mask(regions))
    utils._write_bed_file(args.out_file, resolved, chrom_num)

    return 


if __name__ == '__main__':
    main()
