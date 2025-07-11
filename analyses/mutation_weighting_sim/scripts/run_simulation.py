
import demes
import numpy as np
import pandas 
import msprime
import dpluspy
import sys


chrom = int(sys.argv[1])


# Define dynamic filepaths
rec_map_file = f"sexavg_chr{chrom}.txt.gz"
mut_map_file = f"roulette_1kb_chr{chrom}.csv.gz"
out_fname = f"simulated_chr{chrom}.vcf"


# Define static filepaths
graph_fname = "model.yaml"


# Load a windowed mutation map 
mut_df = pandas.read_csv(mut_map_file)
edges = np.append(
    np.array(mut_df["chromStart"]), np.array(mut_df["chromEnd"])[-1])
vals = np.array(mut_df["mut_map"])
vals[np.isnan(vals)] = 0
mut_map = msprime.RateMap(position=edges, rate=vals)


# Load an empirical recombination map
L = edges[-1]
rec_map = msprime.RateMap.read_hapmap(
    rec_map_file, 
    sequence_length=L, 
    position_col=1, 
    rate_col=2
)


# Load the demography and define the sample configuration
graph = demes.load(graph_fname)
demog = msprime.Demography.from_demes(graph)
samples = {
    "Denisova": 1,
    "Vindija": 1,
    "Yoruba1": 1
}


ts = msprime.sim_ancestry(
    demography=demog,
    samples=samples,
    sequence_length=L,
    recombination_rate=rec_map,
    record_provenance=False,
    discrete_genome=True
)
mts = msprime.sim_mutations(ts, rate=mut_map, record_provenance=False)
with open(out_fname, "w") as fout:
    mts.write_vcf(fout, position_transform=dpluspy.utils._increment1,
        individual_names=["Denisova", "Vindija", "Yoruba1"])


