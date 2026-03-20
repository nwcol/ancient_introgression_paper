"""
Scale a mutation map represented in .csv format to a new average.
"""

import argparse
import numpy as np
import pandas

from dpluspy import utils


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True,
        help="Pathname of input .csv file")
    parser.add_argument("-scale", "--scale", type=float, required=True,
        help="Coefficient by which to scale rates")
    parser.add_argument("-o", "--outfile", required=True, 
        help="Pathname of output .bedgraph file")
    return parser.parse_args()


def scale_windowed_map(infile, scale, outfile):
    input_df = pandas.read_csv(infile)
    mut_map = np.array(input_df["mut_map"], np.float64)
    mut_map *= scale
    data = {
        "chrom": input_df["chrom"],
        "chromStart": input_df["chromStart"], 
        "chromEnd": input_df["chromEnd"],
        "num_sites": input_df["num_sites"],
        "mut_map": mut_map
    }
    pandas.DataFrame(data).to_csv(outfile, index=False, na_rep="nan")
    return


def main():
    args = get_args()
    scale_windowed_map(args.infile, args.scale, args.outfile)
    return


if __name__ == "__main__":
    main()
