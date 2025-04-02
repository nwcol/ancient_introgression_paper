## bootstrap data and save the result in a pickled dictionary

import argparse
import numpy as np
import pickle

from dpluspy import bootstrapping


def get_args():

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--in_files", nargs="*", required=True
    )
    parser.add_argument(
        "-o", "--out_file", required=True
    )
    parser.add_argument(
        "-d", "--denom_files", nargs="*", default=None
    )
    parser.add_argument(
        "--num_reps", type=int, default=1000
    )
    parser.add_argument(
        "--weighted", action="store_true"
    )
    return parser.parse_args()


def main():

    args = get_args()
    if args.denom_files is not None:
        denoms = dict()
        for file in args.denom_files:
            raw = pickle.load(open(file, "rb"))
            for region in raw:
                if region in denoms:
                    raise ValueError(f"{region} occurs twice in denoms")
                denoms[region] = raw[region]
    else:
        denoms = None

    regions = dict()
    for file in args.in_files:
        raw = pickle.load(open(file, "rb"))
        for region in raw:
            if region in regions:
                raise ValueError(f"{region} occurs twice in input")
            if denoms is not None:
                _stats = raw[region]
                if region not in denoms:
                    raise ValueError(f"{region} lacks a denominator")
                denom = denoms[region]
                stats = bootstrapping.exchange_denominators(_stats, denom)
            else:
                stats = raw[region]
            regions[region] = stats

    means, varcovs = bootstrapping.bootstrap(
        regions, 
        num_reps=args.num_reps, 
        # weighted=args.weighted,
        get_reps=False
    )
    bootstrap_data = {}
    bootstrap_data["means"] = means
    bootstrap_data["varcovs"] = varcovs
    # bootstrap_data["reps"] = reps
    bootstrap_data["pop_ids"] = regions[next(iter(regions))]["pop_ids"]
    bootstrap_data["bins"] = regions[next(iter(regions))]["bins"]

    with open(args.out_file, "wb") as fout: pickle.dump(bootstrap_data, fout)

    return 


if __name__ == "__main__":
    main()

