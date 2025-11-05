
import demes
import msprime 
import numpy as np
import dpluspy
import pandas
import sys
import pickle


# Parameters
u_bar = 1.3e-8
pop_ids = ["Yor", "Nea", "Den"]


model = sys.argv[1]
out_fname = sys.argv[2]


def simulate_chrom(chrom):
    rec_map_fname = f"data/recombination_maps/Bherer_b37/sexavg_chr{chrom}.txt.gz"
    mut_map_fname = f"data/mutation_maps/scaled_roulette_b37_10kb/roulette_10kb_chr{chrom}.csv.gz"
    mask_fname = f"data/bed_files/filterbed_1e-4M_exon_buffer/mask_chr{chrom}.bed.gz"
    window_fname = f"data/intervals/1e-4_exon_buffer_1.3Mb/intervals_chr{chrom}.txt"

    # Load mutation map (a .csv file with uniform windows)
    mut_df = pandas.read_csv(mut_map_fname)
    edges = np.append(
        np.array(mut_df["chromStart"]), np.array(mut_df["chromEnd"])[-1])
    vals = np.array(mut_df["mut_map"])
    vals[np.isnan(vals)] = 0
    mut_map = msprime.RateMap(position=edges, rate=vals)

    # Load recombination map
    sequence_length = edges[-1]
    rec_map = msprime.RateMap.read_hapmap(
        rec_map_fname, sequence_length=sequence_length, position_col=1, rate_col=2)

    # Define bins in `r`
    bins = np.logspace(-6, -2, 17)

    # Define demography
    samples = {pop_id: 1 for pop_id in pop_ids}
    graph = demes.load(model)
    demog = msprime.Demography.from_demes(graph)

    # Run the coalescent simulation
    ts = msprime.sim_ancestry(
        demography=demog,
        samples=samples,
        sequence_length=sequence_length,
        recombination_rate=rec_map,
        record_provenance=False,
        discrete_genome=True,
        model=[
            msprime.DiscreteTimeWrightFisher(duration=1000),
            msprime.StandardCoalescent(),
        ]
    )

    # Simulate mutations
    mts = msprime.sim_mutations(ts, rate=mut_map, record_provenance=False)

    # Compute statistics
    stats = dpluspy.parsing.parse_stats(
        mts,
        bed_file=mask_fname,
        rec_map_file=rec_map_fname,
        pos_col="pos",
        map_col="cM",
        mut_map_file=mut_map_fname,
        mut_col="mut_map",
        interval_file=window_fname,
        r_bins=bins,
        chrom=chrom,
        u_bar=u_bar
    )
    return stats


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


sums = {}
for chrom in range(1, 23):
    sums.update(simulate_chrom(chrom))
run_bootstrap(sums, out_fname)

