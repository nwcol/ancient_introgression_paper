# run genome-scale coalescent simulations with or without local mutation
# rate variation.
import demes
import msprime 
import numpy as np
import dpluspy
import pandas
import sys
import pickle


# arguments
apply_lmr_var = bool(int(sys.argv[1]))
out_fstem = sys.argv[2]

print(apply_lmr_var)

# parameters
chroms = list(range(1, 23))
bins = np.concatenate([[0], np.logspace(-6, -1, 21), [0.5]])
u_bar = 1.3e-08
pop_ids = [
    "Altai",
    "Vindija",
    "Chagyrskaya", 
    "Denisova", 
    "Stuttgart",
    "Loschbour",
    "UstIshim",
    "Yoruba1"
]
model_fname = "bherer_full.yaml"
rec_map_fname = lambda n: f"data/recombination_maps/Bherer/sexavg_chr{n}.txt.gz"
mut_map_fname = lambda n: f"data/windowed_mutation_maps/scaled_roulette_1kb/map_chr{n}.csv.gz"
mask_fname = lambda n: f"data/bed_files/1e-4M_buffer/mask_chr{n}.bed.gz"
window_fname = lambda n: f"data/windows/1e-4M_buffer/windows_chr{n}.txt"


def simulate_chrom(chrom):
    # load mutation map (a .csv file with uniform windows)
    mut_df = pandas.read_csv(mut_map_fname(chrom))
    edges = np.append(
        np.array(mut_df["chromStart"]), np.array(mut_df["chromEnd"])[-1])
    vals = np.array(mut_df["mut_map"])
    if apply_lmr_var:
        vals[np.isnan(vals)] = 0
        mut_map = msprime.RateMap(position=edges, rate=vals)
    else:
        # if not applying LMR variation; use the avg rate across the chromosome
        # (number of sites that pass the mask applied)
        num_sites = np.array(mut_df["num_sites"])
        chrom_u = np.nansum(vals * num_sites) / np.sum(num_sites)

    # load recombination map
    sequence_length = edges[-1]
    rec_map = msprime.RateMap.read_hapmap(
        rec_map_fname(chrom), 
        sequence_length=sequence_length, 
        position_col=1, 
        rate_col=2
    )

    # define demography
    samples = {pop_id: 1 for pop_id in pop_ids}
    graph = demes.load(model_fname)
    demog = msprime.Demography.from_demes(graph)

    ts = msprime.sim_ancestry(
        demography=demog,
        samples=samples,
        sequence_length=sequence_length,
        recombination_rate=rec_map,
        record_provenance=False,
        discrete_genome=True)

    # simulate mutations under the proper mutation model
    if apply_lmr_var:
        mts = msprime.sim_mutations(ts, rate=mut_map, record_provenance=False)
    else:
        mts = msprime.sim_mutations(ts, rate=chrom_u, record_provenance=False)

    # compute statistics straight from the tree sequence
    naive_sums = dpluspy.parsing.parse_stats(
        mts,
        bed_file=mask_fname(chrom),
        rec_map_file=rec_map_fname(chrom),
        pos_col="pos",
        map_col="cM",
        interval_file=window_fname(chrom),
        r_bins=bins,
        chrom=chrom
    )
    weighted_sums = dpluspy.parsing.parse_stats(
        mts,
        bed_file=mask_fname(chrom),
        rec_map_file=rec_map_fname(chrom),
        pos_col="pos",
        map_col="cM",
        mut_map_file=mut_map_fname(chrom),
        mut_col="mut_map",
        interval_file=window_fname(chrom),
        r_bins=bins,
        chrom=chrom,
        u_bar=u_bar,
        get_denoms=False
    )
    for key in naive_sums:
        weighted_sums[key]["denoms"] = naive_sums[key]["denoms"]
    return naive_sums, weighted_sums


def run_bootstrap(sums, out_fname):
    bins = sums[next(iter(sums))]["bins"]
    bootreps = dpluspy.bootstrapping.get_bootstrap_reps(sums)
    varcovs = dpluspy.bootstrapping.compute_varcovs(bootreps)
    means = dpluspy.bootstrapping.means_across_regions(sums)
    data = {
        "means": means, 
        "varcovs": varcovs, 
        "bins": bins, 
        "pop_ids": pop_ids
    }
    with open(out_fname, "wb") as fout:
        pickle.dump(data, fout)
    return 


all_naive_sums = {}
all_weighted_sums = {}
for chrom in chroms:
    naive_sums, weighted_sums = simulate_chrom(chrom)
    all_naive_sums.update(naive_sums)
    all_weighted_sums.update(weighted_sums)

with open(f"{out_fstem}_naive_sums.pkl", "wb") as fout:
    pickle.dump(all_naive_sums, fout)
with open(f"{out_fstem}_weighted_sums.pkl", "wb") as fout:
    pickle.dump(all_weighted_sums, fout)

run_bootstrap(all_naive_sums, f"{out_fstem}_naive_stats.pkl")
run_bootstrap(all_weighted_sums, f"{out_fstem}_weighted_stats.pkl")
