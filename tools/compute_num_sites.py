"""
Compute the number of sites in each window under a given window scheme.
"""

import numpy as np
import dpluspy
import sys


def main():
    print("chrom\twindow\tloci\toverlap loci")
    chrom = sys.argv[1]
    window_fname = f"../data/windows/windows.chr{chrom}.txt"
    windows = np.loadtxt(window_fname)
    bed_fname = \
        f"../data/bed_files/filterbed_exons_1e-4M/mask_chr{chrom}.bed.gz"
    positions = dpluspy.utils._read_bed_file_positions(bed_fname) + 1
    for ii, window in enumerate(windows):
        w0, w1, w2 = window
        count0 = np.count_nonzero((positions >= w0) & (positions < w1))
        count1 = np.count_nonzero((positions >= w1) & (positions < w2))
        print(f"{chrom}\t{ii}\t{count0}\t{count1}")
    return


if __name__ == "__main__":
    main()
