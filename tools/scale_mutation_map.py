"""
Scale a site-resolution mutation map by a factor of 1e8, convert it to 16 bit
float entries, and save the result. Should reduce disk space occupied by these
files 4-fold.

Usage: $ python scale_mutation_map.py IN_FNAME OUT_FNAME
"""

import numpy as np
import sys


def scale_mutation_map(in_fname, out_fname):
    mut_map = np.load(in_fname)
    scaled_map = (mut_map * 1e8).astype(np.float16)
    np.save(out_fname, scaled_map)
    return


def main():
    in_fname = sys.argv[1]
    out_fname = sys.argv[2]
    scale_mutation_map(in_fname, out_fname)
    return

if __name__ == "__main__":
    main()
