## make a .bed file specifying sites that are covered in all input bed files

import argparse
import numpy as np

from h2py import utils


def get_args():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--in_file', required=True
    )
    parser.add_argument(
        '-s', '--slop', required=True, type=int
    )
    parser.add_argument(
        '-o', '--out_file', required=True
    )
    return parser.parse_args()


def main():
    #
    args = get_args()
    slop = args.slop
    
    in_regions = utils.read_bedfile(args.in_file)
    chrom_num =  str(np.loadtxt(args.in_file, usecols=0, skiprows=1, dtype=str)[0])

    mask = np.ones(in_regions[-1, 1] + slop, dtype=bool)

    for (start, end) in in_regions:
        mask[start - slop:end + slop] = False
        
    regions = utils.mask_to_regions(mask)
    utils.write_bedfile(args.out_file, regions, chrom_num=chrom_num)

    return 


if __name__ == '__main__':
    main()
