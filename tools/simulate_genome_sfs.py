
import demes
import h2py
from h2py.utils import timestamp
import moments
import msprime
import numpy as np
import pandas
import pickle
import sys


# Input arguments
graph_file = sys.argv[1]
pop_file = sys.argv[2]
out_prefix = sys.argv[3]


# Dynamic inputs.
pop_spec = h2py.matrices._read_pop_file(pop_file)
samples = {pop: 1 for pop in pop_spec}
graph = demes.load(graph_file)
demog = msprime.Demography.from_demes(graph)


# Fixed parameters.
# window size.
ws = int(5e6)
r_bins = np.logspace(-6, -2, 17)
chroms = list(range(1, 23))
# avg mutation rate according to input maps
avg_u = 1.15122e-08
# desired mean u
target_u = 1.3e-8
u_fac = target_u / avg_u
print(timestamp(), f"mutation rate adjustment factor: {u_fac:.6}")


# Filename patterns
mut_map_file = lambda c: f"data/map_chr{c}.tsv.gz"
bed_file = lambda c: f"data/mask_chr{c}.bed.gz"
rec_map_file = lambda c: f"data/sexavg_chr{c}.txt.gz"
genome_file = "data/hg19.genome"
vcf_file = lambda c: f"data_chr{c}.vcf"


L_df = pandas.read_csv(genome_file, names=["chrom", "L"], sep="\\s+")
L_map = {int(c.replace("chr", "")): next(iter(L_df[L_df["chrom"] == c]["L"]))
         for c in L_df["chrom"]}

interval_map = {c: [[i*ws, (i+1)*ws] for i in range(int(L_map[c]//ws) + 1)]
                for c in chroms}


def simulate_chromosome(chrom):
    """``chrom`` must be an integer."""
    # Setup
    graph = demes.load(graph_file)
    L = L_map[chrom]

    rec_map = msprime.RateMap.read_hapmap(rec_map_file(chrom), sequence_length=L)

    mut_tbl = pandas.read_table(mut_map_file(chrom), engine="python", sep=None)
    mut_rate = np.array(mut_tbl["mutRate"]) * u_fac
    mut_rate[np.isnan(mut_rate)] = 0
    avg_rate = np.mean(mut_rate[mut_rate > 0])
    print(timestamp(), f"average mutation rate on {chrom}: {avg_rate:.6}")
    starts = np.array(mut_tbl["chromStart"])
    mut_pos = np.append(starts[starts < L], L)
    mut_map = msprime.RateMap(position=mut_pos, rate=mut_rate[starts < L])

    ts = msprime.sim_ancestry(
        samples,
        demography=demog,
        sequence_length=L,
        recombination_rate=rec_map,
        record_provenance=False
    )
    ts = msprime.sim_mutations(ts, rate=mut_map, record_provenance=False)
    with open(vcf_file(chrom), "w+") as f:
        ts.write_vcf(
            f, contig_id=chrom, position_transform=lambda x: 1 + np.array(x)),
            individual_names=list(samples.keys()))
    return


def compute_stats(chrom):
    """calculate H2 sums for a chromosome"""
    chrom_sums = {}
    chrom_fs = {}
    intervals = interval_map[chrom]
    for i, interval in enumerate(intervals):
        interval_sums = h2py.parsing.compute_h2_stats(
            vcf_file=vcf_file(chrom),
            pop_file=pop_file,
            interval=interval,
            bed_file=bed_file(chrom),
            rec_map_file=rec_map_file(chrom),
            rec_pos_col="pos",
            rec_map_col="cM",
            mut_map_file=mut_map_file(chrom),
            r_bins=r_bins,
            report=True,
        )
        if interval_sums["denoms"][-1] == 0:
            print(timestamp(), f"chromosome {chrom} interval {interval} has 0 sites: passing")
            continue
        chrom_sums[f"{chrom}_{i}"] = interval_sums
        # An error will be raised by moments if the window has no polymorphism
        # (e.g., if the number of accessible sites is very low). Avoid this:
        if np.sum(interval_sums["sums"][-1]) == 0:
            print(timestamp(), f"chromosome {chrom} interval {interval} has no polymorphism: passing")
            continue
        # Ref is interpreted as ancestral (this is correct, as msprime uses
        # this convention).
        chrom_fs[f"{chrom}_{i}"] = moments.Spectrum.from_vcf(
            vcf_file(chrom),
            pop_file=pop_file,
            bed_file=bed_file(chrom),
            interval=interval,
        )
        print(timestamp(), f"parsed stats on chromosome {c} interval {interval}")

    return chrom_sums, chrom_fs


print(timestamp(), "running simulations")

all_h2 = {}
all_fs = {}

for c in chroms:
    simulate_chromosome(c)
    print(timestamp(), f"simulated chromosome {c}")
    chrom_h2, chrom_fs = compute_stats(c)
    all_h2.update(chrom_h2)
    all_fs.update(chrom_fs)

# Boostrap H2 statistics.
u_bar, n_sites = h2py.parsing.compute_avg_u(all_h2)
u_bar *= u_fac  # adjust the measured avg mutation rate
print(timestamp(), f"u_bar = {u_bar:.6}")
print(timestamp(), f"L = {n_sites}")
fac = u_bar ** 2
for label in all_h2:
    for i in range(len(all_h2[label]["sums"]) - 1):
        all_h2[label]["sums"][i] *= fac

boot_h2 = h2py.parsing.bootstrap_data(all_h2)

with open(f"{out_prefix}_boot_h2.pkl", "wb") as fout:
    pickle.dump(boot_h2, fout)

# Process SFS statistics.
template = all_fs[next(iter(all_fs))]
sum_fs = np.sum([all_fs[label] for label in all_fs], axis=0)
sum_fs = moments.Spectrum(sum_fs, pop_ids=template.pop_ids)
fs_out = {"fs": sum_fs, "L": n_sites}

with open(f"{out_prefix}_sum_fs.pkl", "wb") as fout:
    pickle.dump(fs_out, fout)

# Save interval-specific raw outputs.
#with open(f"{out_prefix}_all_h2.pkl", "wb") as fout:
#    pickle.dump(all_h2, fout)

#with open(f"{out_prefix}_all_fs.pkl", "wb") as fout:
#    pickle.dump(all_fs, fout)

print(timestamp(), "genome simulation complete!")

