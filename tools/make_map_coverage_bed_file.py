"""
Create a BED file representing non-NaN coverage of a site-resolution
mutation map (.npy format).
"""

import argparse
import gzip
import numpy as np

import dpluspy 


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True,
        help='Pathname of input .npy file')
    parser.add_argument('-chrom_num', '--chrom_num', required=True,
        help='String representation of the chromosome number')
    parser.add_argument('-o', '--outfile', required=True,
        help='Pathname of output .bed file')
    return parser.parse_args()


def make_map_coverage_mask(infile, outfile, chrom_num):
    mut_map = np.load(infile)
    bool_mask = np.isnan(mut_map)
    regions = dpluspy.utils._mask_to_regions(bool_mask)
    dpluspy.utils._write_bed_file(outfile, regions, chrom_num)
    return 


def main():
    args = get_args()
    make_map_coverage_mask(args.infile, args.outfile, args.chrom_num)
    return


if __name__ == "__main__":
    main()

