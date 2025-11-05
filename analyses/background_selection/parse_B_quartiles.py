
import pickle
import sys
import numpy as np
from dpluspy import parsing


# ARGUMENTS
chrom = int(sys.argv[1])
quartile = sys.argv[2]


out_file = f"stats_B_strat_quartile_{quartile}_chr{chrom}.pkl"


# Define dynamic filepaths
vcf_file = f"merged.variants.chr{chrom}.vcf.gz"
bed_file = f"mask_chr{chrom}.bed.gz"
rec_map_file = f"sexavg_chr{chrom}.txt.gz"
mut_map_file = f"roulette_b37_chr{chrom}.npy"
interval_file = f"arms_chr{chrom}.txt"
pop_file = "populations_sub"

# mutation normalization factor
u_bar = 1.1221531631721093e-08


# Define bins
bins = np.concatenate(([0], np.logspace(-7, -1, 25), 
    [1.77827941e-01, 3.16227766e-01, 4.99999e-01]))


stats = parsing.parse_stats(
    vcf_file,
    bed_file=bed_file,
    pop_file=pop_file,
    rec_map_file=rec_map_file,
    pos_col="pos",
    map_col="cM",
    u_bar=u_bar,
    r_bins=bins,
    mut_map_file=mut_map_file,
    interval_file=interval_file,
    chrom=chrom
)


with open(out_file, "wb") as fout:
    pickle.dump(stats, fout)
