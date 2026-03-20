"""
Intersect several BED files and write a single output BED file recording the
intervals where all inputs had coverage.
"""

import argparse
import numpy as np

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infiles', nargs='*', required=True,
        help='Input BED filepaths')
    parser.add_argument('-o', '--outfile', required=True,
        help='Output BED filepath')
    return parser.parse_args()


def intersect_bed_files(infiles, outfile):
    """
    Load BED files, form their intersection, and save it.
    """
    region_arrs = []
    chrom_nums = []
    for file in infiles:
        regions, chrom = utils._read_bed_file(file)
        region_arrs.append(regions)
        chrom_nums.append(chrom)
    chrom_num = chrom_nums[0]
    intersection = utils._intersect_regions(region_arrs)
    utils._write_bed_file(outfile, intersection, chrom_num)
    return


def main():
    args = get_args()
    intersect_bed_files(args.infiles, args.outfile)
    return 


if __name__ == '__main__':
    main()
