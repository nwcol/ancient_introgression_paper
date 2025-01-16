## make a .bed file specifying sites that are covered in all input bed files

import argparse
import numpy as np

from h2py import util


def get_args():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--in_file', required=True
    )
    parser.add_argument(
        '-s', '--subtrahend_file', required=True
    )
    parser.add_argument(
        '-o', '--out_file', required=True
    )
    return parser.parse_args()


def main():
    #
    args = get_args()
    
    in_regions = util.read_bedfile(args.in_file)
    sub_regions = util.read_bedfile(args.subtrahend_file)
    chrom_num =  str(np.loadtxt(args.in_file, usecols=0, skiprows=1, dtype=str)[0])

    max_length = max([in_regions[-1, 1], sub_regions[-1, 1]])
    # zeros are `in` the mask (included); ones are outside it
    mask = np.ones(max_length, dtype=bool)

    for (start, end) in in_regions:
        mask[start:end] = False

    for (start, end) in sub_regions:
        mask[start:end] = True  
        
    regions = util.mask_to_regions(mask)
    util.write_bedfile(args.out_file, regions, chrom_num=chrom_num)

    return 


if __name__ == '__main__':
    main()
