
import dpluspy
import pickle
import numpy as np
import sys


out_fname = sys.argv[1]
in_fpattern = sys.argv[2]


def run_bootstrap(in_fpattern, out_fname):

    denom_fnames = [f"../../one_genome/raw_stats/prelim_chr_{c}_stats2.pkl" 
                    for c in range(1, 23)]
    denoms = {}
    for fname in denom_fnames:
        with open(fname, "rb") as fin:
            data = pickle.load(fin)
            for (x, y) in data:
                denoms[(int(x), int(y))] = data[(x, y)]

    in_fnames = [f"{in_fpattern}_chr_{c}_stats.pkl" for c in range(1, 23)]
    sums = {}
    for fname in in_fnames:
        with open(fname, "rb") as fin:
            data = pickle.load(fin)
            for (x, y) in data:
                sums[(int(x), int(y))] = data[(x, y)]

    for region in sums:
        sums[region]["denoms"] = denoms[region]["denoms"]

    bins = sums[next(iter(sums))]["bins"]

    # I misspecified these and need to manually replace them
    pop_ids = ["Nea", "Yor"]

    replicates = dpluspy.bootstrapping.get_bootstrap_reps(sums, weighted=False)
    varcovs = dpluspy.bootstrapping.compute_varcovs(replicates)
    means = dpluspy.bootstrapping.means_across_regions(sums)

    data = {
        "means": means, 
        "varcovs": varcovs, 
        "bins": bins, 
        "pop_ids": pop_ids,
        "replicates": replicates
    }
    with open(out_fname, "wb") as fout:
        pickle.dump(data, fout)

    return 


def main():
    run_bootstrap(in_fpattern, out_fname)
    return


main()

