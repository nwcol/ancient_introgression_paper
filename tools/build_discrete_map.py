"""
Take a site-resolution mutation map (.npy format) and build a windowed analog.
"""

import argparse
import numpy as np

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--in_file', required=True,
        help='Pathnames of input .npy file')
    parser.add_argument('-scale', '--scale', type=int,
        help='Window scale, in bp')
    parser.add_argument('-chrom', '--chrom', default='0', 
        help='Chromosome number for output file')
    parser.add_argument('-o', '--out_file', required=True, 
        help='Pathname of output .bedgraph file')
    return parser.parse_args()


def main():
    args = get_args()
    mut_map = np.load(args.in_file)
    L = len(mut_map)
    scale = args.scale
    windows = np.stack(
        (np.arange(0, L, scale), np.arange(scale, L + scale, scale)), axis=1
    )
    window_map = np.zeros(len(windows), dtype=np.float64)
    fill_val = np.nanmean(mut_map)
    for i, (start, end) in enumerate(windows):
        if np.all(np.isnan(mut_map[start:end])):
            window_map[i] = fill_val 
        else:
            window_map[i] = np.nanmean(mut_map[start:end])
    data = {'map': window_map}
    utils._write_bedgraph_file(args.out_file, windows, data, args.chrom)

    return


if __name__ == '__main__':
    main()
