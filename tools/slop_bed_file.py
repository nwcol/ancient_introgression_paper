"""
Extend the regions in a BED file by some physical distance and save the result.
"""

import argparse
import h2py


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, 
        help='Pathname of input BED file')
    parser.add_argument('-slop', '--slop', required=True, type=int,
        help='Slop distance, in bp')
    parser.add_argument('-o', '--outfile', required=True,
        help='Pathname of output BED file')
    return parser.parse_args()


def slop_bed_file(infile, slop, outfile):
    mask = h2py.GeneticMask.from_bed_file(infile)
    mask = mask.add_flank(slop)
    mask.to_bed_file(outfile)
    return


def main():
    args = get_args()
    slop_bed_file(args.infile, args.slop, args.outfile)
    return


if __name__ == '__main__':
    main()
