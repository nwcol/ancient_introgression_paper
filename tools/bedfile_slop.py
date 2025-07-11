"""
Extend the regions in a BED file by some physical distance and save the result.
"""

import argparse

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', required=True, 
        help='Pathname of input BED file')
    parser.add_argument('-buffer', '--buffer', required=True, type=int,
        help='Slop distance, in bp')
    parser.add_argument('-o', '--out_file', required=True,
        help='Pathname of output BED file')
    return parser.parse_args()


def main():
    args = get_args()
    regions, chrom_num = utils._read_bed_file(args.in_file)
    buffer = args.buffer
    regions[:, 0] -= buffer
    regions[:, 1] += buffer
    regions[regions < 0] = 0
    resolved = utils._mask_to_regions(utils._regions_to_mask(regions))
    utils._write_bed_file(args.out_file, resolved, chrom_num)
    return 


if __name__ == '__main__':
    main()
