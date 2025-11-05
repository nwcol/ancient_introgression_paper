
import pickle
import sys
import numpy as np
from dpluspy import parsing


# Arguments
chrom = int(sys.argv[1])
map_name = str(sys.argv[2])
rec_map_file = str(sys.argv[3])
pos_col = str(sys.argv[4])
map_col = str(sys.argv[5])


out_file = f"sums_{map_name}_chr{chrom}.pkl"


# Define dynamic filepaths
vcf_file = f"merged.variants.chr{chrom}.vcf.gz"
bed_file = f"mask_chr{chrom}.bed.gz"
mut_map_file = f"roulette_b37_chr{chrom}.npy"
interval_file = f"arms_chr{chrom}.txt"
vcf_file = f"merged.variants.chr{chrom}.vcf.gz"
pop_file = "populations"

# mutation normalization factor
u_bar = 1.1221531631721093


# Define bins
bins = np.concatenate(([0], np.logspace(-7, -1, 25), 
    [1.77827941e-01, 3.16227766e-01, 4.99999e-01]))


stats = parsing.parse_stats(
    vcf_file,
    bed_file=bed_file,
    label_file=pop_file,
    rec_map_file=rec_map_file,
    pos_col=pos_col,
    map_col=map_col,
    u_bar=u_bar,
    r_bins=bins,
    mut_map_file=mut_map_file,
    interval_file=interval_file,
    chrom=chrom
)


with open(out_file, "wb") as fout:
    pickle.dump(stats, fout)

