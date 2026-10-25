"""
Create a BED file representing non-NaN coverage of a site-resolution
mutation map (.npy format).
"""

import argparse
import h2py
import gzip
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True,
        help='Pathname of input .npy file')
    parser.add_argument('-chrom_num', '--chrom_num', required=True,
        help='String representation of the chromosome number')
    parser.add_argument('-o', '--outfile', required=True,
        help='Pathname of output .bed file')
    return parser.parse_args()


def make_map_coverage_mask(infile, outfile, chrom):
    mut_map = np.load(infile)
    bool_mask = np.logical_not(np.isnan(mut_map))
    mask = h2py.GeneticMask(bool_mask, chrom=chrom)
    mask.to_bed_file(outfile)
    return


def main():
    args = get_args()
    make_map_coverage_mask(args.infile, args.outfile, args.chrom_num)
    return


if __name__ == "__main__":
    main()

