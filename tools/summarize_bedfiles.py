
import argparse
import numpy as np
import sys

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_files', nargs='*', default=[])
    parser.add_argument('-check', '--check', action='store_true',
        help='If used, removes any overlap between regions when counting sites')
    return parser.parse_args()


def main():
    args = get_args()
    filenames = args.in_files
    check = bool(args.check)
    tot = 0
    for filename in filenames:
        regions, _ = utils._read_bed_file(filename)
        if check:
            regions = utils._mask_to_regions(utils._regions_to_mask(regions))
        num_sites = int(np.sum(np.diff(regions, axis=1)))
        tot += num_sites
        print(f"{filename}\t{num_sites}")
    print(f"total\t{tot}")

    return


if __name__ == "__main__":
    main()
