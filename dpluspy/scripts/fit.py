## fit D+ to data.

import argparse

from dpluspy import inference


def get_args():
    
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-d", "--data_file", type=str, required=True
    )
    parser.add_argument(
        "-g", "--graph_file", type=str, required=True
    )
    parser.add_argument(
        "-p", "--param_file", type=str, required=True
    )
    parser.add_argument(
        "-u", "--u", type=float, required=True
    )
    parser.add_argument(
        "-o", "--out_file", type=str, required=True
    )

    # optional arguments
    parser.add_argument(
        "-perturb", "--perturb", type=float,  default=0
    )    
    parser.add_argument(
        "-log", "--log", action="store_true",
        help="optimize over the logarithm of the parameters"
    )
    parser.add_argument(
        "-v", "--verbose", type=int, default=1
    )
    parser.add_argument(
        "-m", "--max_iter", type=int, default=1000
    )
    parser.add_argument(
        "--method", type=str, default="fmin"
    )
    parser.add_argument(
        "-one_locus", "--one_locus", action="store_true"
    )

    return parser.parse_args()


def main():

    args = get_args()
    data = inference.load_statistics(args.data_file, args.graph_file)
    pop_ids, bins, means, varcovs = data
    inference.optimize(
        args.graph_file,
        args.param_file,
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=args.u,
        verbose=args.verbose,
        method=args.method,
        max_iter=args.max_iter,
        log=args.log,
        out_file=args.out_file,
        perturb=args.perturb,
        one_locus=args.one_locus
    )

    return


if __name__ == "__main__":
    main()
    