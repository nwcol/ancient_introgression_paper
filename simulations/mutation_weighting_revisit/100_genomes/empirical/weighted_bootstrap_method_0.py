
import dpluspy
import pickle
import numpy as np
import sys


out_fname = sys.argv[1]
in_fnames = sys.argv[2]


def run_bootstrap(in_fpattern, out_fname):

    denom_fnames = [f"../../one_genome/raw_stats/prelim_chr_{c}_stats2.pkl" 
                    for c in range(1, 23)]
    denoms = {}
    for fname in denom_fnames:
        with open(fname, "rb") as fin:
            denoms.update(pickle.load(fin))

    in_fnames = [f"{in_fpattern}_chr_{c}_stats.pkl" for c in range(1, 23)]
    sums = {}
    for fname in in_fnames:
        with open(fname, "rb") as fin:
            sums.update(pickle.load(fin))
    bins = sums[next(iter(sums))]["bins"]
    # I misspecified these; need to replace them
    # pop_ids = sums[next(iter(sums))]["pop_ids"]
    pop_ids = ["Nea", "Yor"]

    for region in sums:
        sums[region]["denoms"] = denoms[region]["denoms"]
        sums[region]["mut_facs"] = denoms[region]["mut_facs"]

    replicates = dpluspy.bootstrapping.get_bootstrap_reps(
        sums, weighted=True, aggregate=False)
    varcovs = dpluspy.bootstrapping.compute_varcovs(replicates)
    means = dpluspy.bootstrapping.weighted_means_across_regions(sums)

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
    run_bootstrap(in_fnames, out_fname)
    return


main()

