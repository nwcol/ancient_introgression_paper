
import numpy as np
import sys

from dpluspy import utils


def main():

    in_file = sys.argv[1]
    out_file = sys.argv[2]
    mut_map = np.load(in_file)
    L = len(mut_map)
    scale = 1000
    windows = np.stack(
        (np.arange(0, L, scale), np.arange(scale, L + scale, scale)), axis=1
    )
    discrete_map = np.zeros(len(windows), dtype=np.float64)
    for i, (start, end) in enumerate(windows):
        if np.all(np.isnan(mut_map[start:end])):
            continue 
        else:
            discrete_map[i] = np.nanmean(mut_map[start:end])
    utils.write_bedgraph(out_file, "chr", windows, {"u": discrete_map})

    return


main()
