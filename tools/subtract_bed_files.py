"""
Take the set difference of the sites in intervals of --in_file1 and in 
--in_file2. 
"""

import argparse
import dpluspy
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--infile1', required=True,
        help='BED file to subtract from (minuend)')
    parser.add_argument('-i2', '--infile2', required=True,
        help='BED file to subtract (subtrahend)')
    parser.add_argument('-o', '--outfile', required=True,
        help='output BED file')
    return parser.parse_args()


def subtract_bed_files(infile1, infile2, outfile):
    """
    Compute the difference as the intersection of sites in the first file
    and the complement of sites in the second file.
    """
    regions1, chrom_num = dpluspy.utils._read_bed_file(infile1)
    regions2, _ = dpluspy.utils._read_bed_file(infile2)
    max_length = max([regions1[-1, 1], regions2[-1, 1]])
    comp_regions2 = dpluspy.utils._mask_to_regions(
        ~dpluspy.utils._regions_to_mask(regions2, length=max_length)
    )
    difference = dpluspy.utils._intersect_regions([regions1, comp_regions2])
    dpluspy.utils._write_bed_file(outfile, difference, chrom_num)
    return 


def main():
    args = get_args()
    subtract_bed_files(args.infile1, args.infile2, args.outfile)
    return


if __name__ == '__main__':
    main()
