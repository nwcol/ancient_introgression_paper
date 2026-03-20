"""
Extend the regions in a BED file by some physical distance and save the result.
"""

import argparse
import dpluspy


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, 
        help='Pathname of input BED file')
    parser.add_argument('-buffer', '--buffer', required=True, type=int,
        help='Slop distance, in bp')
    parser.add_argument('-o', '--outfile', required=True,
        help='Pathname of output BED file')
    return parser.parse_args()


def slop_bed_file(infile, buffer, outfile):
    args = get_args()
    regions, chrom_num = dpluspy.utils._read_bed_file(infile)
    regions[:, 0] -= buffer
    regions[:, 1] += buffer
    regions[regions < 0] = 0
    resolved = dpluspy.utils._mask_to_regions(
        dpluspy.utils._regions_to_mask(regions))
    dpluspy.utils._write_bed_file(outfile, resolved, chrom_num)
    return 


def main():
    args = get_args()
    slop_bed_file(args.infile, args.buffer, args.outfile)
    return


if __name__ == '__main__':
    main()
