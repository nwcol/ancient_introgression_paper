

import dpluspy
import pickle
import numpy as np
import sys


in_fnames = sys.argv[1:]


def compute_u(in_fnames):

    print(f"Running bootstrap on {len(in_fnames)} files")
    sums = {}
    for fname in in_fnames:
        with open(fname, "rb") as fin:
            sums.update(pickle.load(fin))

    numer = np.sum([sums[x]["mut_facs"][-1] for x in sums])
    denom = np.sum([sums[x]["denoms"][-1] for x in sums])
    print(numer, denom)
    print(numer / denom)

    return 


def compute_two_locus_u(in_fnames):

    print(f"Running bootstrap on {len(in_fnames)} files")
    sums = {}
    for fname in in_fnames:
        with open(fname, "rb") as fin:
            sums.update(pickle.load(fin))

    numer = np.sum([sums[x]["mut_facs"][:-1] for x in sums])
    denom = np.sum([sums[x]["denoms"][:-1] for x in sums])
    print(numer, denom)
    print((numer / denom) ** 0.5)

    return 


def main():
    compute_u(in_fnames)
    compute_two_locus_u(in_fnames)
    return


main()

