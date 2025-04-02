## make a .bed file specifying sites that are covered in all input bed files

import argparse
import numpy as np

from h2py import utils


def get_args():
    # get args
    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-i', '--in_files', nargs='*', required=True
    )
    parser.add_argument(
        '-o', '--out_file', required=True
    )
    return parser.parse_args()


def main():
    #
    args = get_args()

    regionss = []
    
    for file in args.in_files:
        regions = utils.read_bedfile(file)
        regionss.append(regions)

    chrom_nums = []
    for file in args.in_files:
        file_chroms = np.loadtxt(file, usecols=0, skiprows=1, dtype=str)
        chrom_nums += list(set(list(file_chroms)))

    chrom_num = list(set(chrom_nums))[0]

    max_length = max([regions[-1, 1] for regions in regionss])
    coverage = np.zeros(max_length, dtype=np.int8) # don't intersect more than 128 masks

    for regions in regionss:
        for (start, end) in regions:
            coverage[start:end] += 1
        
    mask = coverage < len(regionss)
    intersected_regions = utils.mask_to_regions(mask)

    utils.write_bedfile(args.out_file, intersected_regions, chrom_num=chrom_num)
    return 


if __name__ == '__main__':
    main()
