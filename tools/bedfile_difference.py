"""
Take the set difference of the sites in intervals of --in_file1 and in 
--in_file2. 
"""

import argparse
import numpy as np

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--in_file1', required=True,
        help='BED file to s')
    parser.add_argument('-i2', '--in_file2', required=True,
        help='BED file to s')
    parser.add_argument('-out', '--out_file', required=True,
        help='Output BED file')
    return parser.parse_args()


def main():
    """
    Compute the difference as the intersection of sites in the first file
    and the complement of sites in the second file.
    """
    args = get_args()
    regions1, chrom_num = utils._read_bed_file(args.in_file1)
    regions2, _ = utils._read_bed_file(args.in_file2)
    max_length = max([regions1[-1, 1], regions2[-1, 1]])
    comp_regions2 = utils._mask_to_regions(
        ~utils._regions_to_mask(regions2, length=max_length)
    )
    difference = utils._intersect_regions([regions1, comp_regions2])
    utils._write_bed_file(args.out_file, difference, chrom_num)

    return 


if __name__ == '__main__':
    main()
