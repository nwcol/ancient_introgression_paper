"""
Intersect several BED files and write a single output BED file recording the
intervals where all inputs had coverage.
"""

import argparse
import h2py
import numpy as np


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inputs', nargs='*', required=True,
        help='Input BED filepaths')
    parser.add_argument('-o', '--output', required=True,
        help='Output BED filepath')
    return parser.parse_args()


def intersect_bed_files(inputs, output):
    """
    Load BED files, form their intersection, and save it.
    """
    masks = [h2py.GeneticMask.from_bed_file(f) for f in inputs]
    result = h2py.masks.intersect_masks(masks)
    result.to_bed_file(output)
    return


def main():
    args = get_args()
    intersect_bed_files(args.inputs, args.output)
    return 


if __name__ == '__main__':
    main()
