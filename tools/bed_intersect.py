"""
Intersect several BED files and write a single output BED file recording the
intervals where all inputs had coverage.
"""

import argparse
import numpy as np

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_files', nargs='*', required=True,
        help='Input BED files')
    parser.add_argument('-o', '--out_file', required=True,
        help='Output BED file')
    return parser.parse_args()


def main():
    """
    Load BED files, form their intersection, and save it.
    """
    args = get_args()
    region_arrs = []
    chrom_nums = []
    for file in args.in_files:
        regions, chrom = utils._read_bed_file(file)
        region_arrs.append(regions)
        chrom_nums.append(chrom)
    chrom_num = chrom_nums[0]
    intersection = utils._intersect_regions(region_arrs)
    utils._write_bed_file(args.out_file, intersection, chrom_num)

    return 


if __name__ == '__main__':
    main()
