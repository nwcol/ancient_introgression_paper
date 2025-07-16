# Parse D+ on the cluster.

import pickle
import sys
import numpy as np
from dpluspy import parsing


# ARGUMENTS
chrom = int(sys.argv[1])
map_name = str(sys.argv[2])
map_file = str(sys.argv[3])
pos_col = str(sys.argv[4])
map_col = str(sys.argv[5])


out_file = f"{map_name}_stats_chr{chrom}.pkl"


# Define dynamic filepaths
mutation_map_file = f"roulette_1kb_chr{chrom}.csv.gz"
bed_file = f"mask_chr{chrom}.bed.gz"
regions = f"windows.chr{chrom}.txt"
vcf_file = f"Eurasians.variants.chr{chrom}.vcf.gz"
pop_file = "populations"


# Define bins
bins = np.concatenate(([0], np.logspace(-6, -1, 21), [0.1778, 0.3162, 0.4999]))


stats = parsing.parse_statistics(
    vcf_file,
    bed_file,
    pop_file=pop_file,
    rec_map_file=map_file,
    pos_col=pos_col,
    map_col=map_col,
    r_bins=bins,
    mut_map_file=mutation_map_file,
    mut_map_col="mut_map",
    regions_file=regions,
    chrom=chrom
)


with open(out_file, "wb") as fout:
    pickle.dump(stats, fout)
