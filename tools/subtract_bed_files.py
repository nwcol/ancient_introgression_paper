
import argparse
import h2py
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i1', '--input0', required=True,
        help='BED file to subtract from (minuend)')
    parser.add_argument('-i2', '--input1', required=True,
        help='BED file to subtract (subtrahend)')
    parser.add_argument('-o', '--output', required=True,
        help='output BED file')
    return parser.parse_args()


def subtract_bed_files(input0, input1, output):
    mask0 = h2py.GeneticMask.from_bed_file(input0)
    mask1 = h2py.GeneticMask.from_bed_file(input1)
    result = h2py.masks.subtract_masks(mask0, mask1)
    result.to_bed_file(output)
    return


def main():
    args = get_args()
    subtract_bed_files(args.input0, args.input1, args.output)
    return


if __name__ == '__main__':
    main()
