
import pickle
import sys
import numpy as np
from dpluspy import parsing


# ARGUMENTS
chrom = int(sys.argv[1])
interval_file = sys.argv[2]
overhang = bool(sys.argv[3])
out_file = sys.argv[4]


# Define dynamic filepaths
vcf_file = f"simulated_chr{chrom}.vcf"
bed_file = f"mask_chr{chrom}.bed.gz"
rec_map_file = f"sexavg_chr{chrom}.txt.gz"
mut_map_file = f"roulette_1kb_chr{chrom}.csv.gz"


# Define bins
bins = np.concatenate(([0], np.logspace(-7, -1, 25), 
    [1.77827941e-01, 3.16227766e-01, 4.99999e-01]))

pop_mapping = {
    "Denisova": ["Denisova"],
    "Vindija": ["Vindija"],
    "Yoruba1": ["Yoruba1"]
}

stats = parsing.parse_stats(
    vcf_file,
    bed_file=bed_file,
    pop_mapping=pop_mapping,
    rec_map_file=rec_map_file,
    pos_col="pos",
    map_col="cM",
    r_bins=bins,
    mut_map_file=mut_map_file,
    mut_col="mut_map",
    interval_file=interval_file,
    chrom=chrom,
    overhang=overhang
)


with open(out_file, "wb") as fout:
    pickle.dump(stats, fout)
