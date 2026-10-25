"""
Compute estimated observed H2 for a single chromosome.
"""

import argparse
import demes
import h2py
from h2py.utils import timestamp
import moments
import msprime
import numpy as np
import pandas
import pickle


def get_args():
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument(
        "-v",
        "--vcf",
        required=True,
        help="VCF file"
    )
    parser.add_argument(
        "-b",
        "--bed",
        required=True,
        help="BED accessibility mask file"
    )
    parser.add_argument(
        "-p",
        "--pops",
        required=True,
        help="Population file"
    )
    parser.add_argument(
        "-L",
        "--L",
        type=int,
        required=True,
        help="Sequence length"
    )
    parser.add_argument(
        "-r",
        "--rec_map",
        required=True,
        help="Hapmap-format recombination map"
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output path (.pkl)"
    )

    # optional arguments with defaults
    parser.add_argument(
        "--rec_pos_col",
        default="pos",
        help="Recombination map position column"
    )
    parser.add_argument(
        "--rec_map_col",
        default="cM",
        help="Recombination map genetic map column"
    )
    parser.add_argument(
        "--mut_map",
        default=None,
        help="Mutation map"
    )
    parser.add_argument(
        "--window_size",
        type=int,
        default=5000000,
        help="Length for intervals"
    )
    parser.add_argument(
        "--report",
        action="store_true",
    )
    parser.add_argument(
        "--bins_start",
        type=float,
        default=1e-6
    )
    parser.add_argument(
        "--bins_end",
        type=float,
        default=1e-2
    )
    parser.add_argument(
        "--n_bins",
        type=int,
        default=17
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    r_bins = np.logspace(
        np.log10(args.bins_start), np.log10(args.bins_end), args.n_bins)
    ws = args.window_size
    intervals = [[i * ws, (i + 1) * ws] for i in range(int(args.L//ws) + 1)]
    sums = {}
    for i, interval in enumerate(intervals):
        interval_sums = h2py.parsing.compute_h2_stats(
            vcf_file=args.vcf,
            pop_file=args.pops,
            interval=interval,
            bed_file=args.bed,
            rec_map_file=args.rec_map,
            rec_pos_col=args.rec_pos_col,
            rec_map_col=args.rec_map_col,
            mut_map_file=args.mut_map,
            r_bins=r_bins,
            report=True,
        )
        if interval_sums["denoms"][-1] == 0:
            print(timestamp(), f"{args.vcf} interval {interval} has 0 sites: passing")
            continue
        sums[f"{args.vcf}_{i}"] = interval_sums
        print(timestamp(), f"parsed stats from {args.vcf} interval {interval}")

    with open(args.out, "wb") as fout:
        pickle.dump(sums, fout)
    print(timestamp(), f"wrote stats from {len(sums)} intervals to {args.out}")
