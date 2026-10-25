"""
Fit a model to observed H2 statistics.
"""

import argparse
import demes
import h2py
from h2py.utils import timestamp
import sys
import time
import numpy as np
import pickle


def parse_args():
    parser = argparse.ArgumentParser()
    # required arguments
    parser.add_argument(
        "-g",
        "--graph",
        required=True,
        help="Demes-format YAML file with initial parameter guesses"
    )
    parser.add_argument(
        "-p",
        "--params",
        required=True,
        help="Moments-Demes style YAML file specifying parameters"
    )
    parser.add_argument(
        "-d",
        "--data",
        required=True,
        help="Input data file (.pkl; see module docstring for details)"
    )
    parser.add_argument(
        "-o",
        "--out",
        required=True,
        help="Output path for fitted model"
    )

    # optional arguments, with default values
    parser.add_argument(
        "--sample",
        type=int,
        default=None,
        help="Bootstrap replicate to use (if `data` is a set of bootstramp samples)"
    )
    parser.add_argument(
        "-u",
        "--u",
        type=float,
        default=1.3e-8,
        help="Mutation rate (default 1.3e-8)"
    )
    parser.add_argument(
        "--perturb",
        type=float,
        default=0,
        help="Random perturbation to initial parameter values (default 0)"
    )
    parser.add_argument(
        "--jitter",
        type=float,
        default=0.01,
        help="Random perturbation between rounds of fitting (default 0.01)"
    )
    parser.add_argument(
        "--report",
        type=int,
        default=10,
        help="Interval for status reports (default 10). 0 suppresses verbosity"
    )
    parser.add_argument(
        "--max_iter",
        type=int,
        default=1000,
        help="Max iterations for each fitting stage (default 1000)"
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1e-2,
        help="Absolute threshold in LL for cessation of fitting"
    )
    parser.add_argument(
        "--max_runtime",
        type=float,
        default=None,
        help="Limit on runtime (s)."
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    with open(args.data, "rb") as fin:
        data = pickle.load(fin)

    if args.sample is not None and args.sample >= 0:
        data = {
            "means": data["bootreps"][args.sample],
            "varcovs": data["varcovs"],
            "bins": data["bins"],
            "pops": data["pops"],
        }

    graph = demes.load(args.graph)
    deme_names = [d.name for d in graph.demes]
    pops = [p for p in data["pops"] if p in deme_names]
    data = h2py.parsing.subset_data(data, to_pops=pops)

    t0 = time.time()
    converged = False

    kwargs = {
        "pops": data["pops"],
        "r_bins": data["bins"],
        "u": args.u,
        "use_log": True,
        "max_iter": args.max_iter,
        "report": args.report,
        "output": args.out,
        "overwrite": True
    }

    print(f"{timestamp()} BEGINNING INITIAL FMIN PHASE")
    _, __, ll0 = h2py.inference.optimize(
        args.graph,
        args.params,
        data["means"],
        data["varcovs"],
        method="fmin",
        perturb=args.perturb,
        **kwargs
    )

    input = args.out
    converged = False
    counter = 0
    ll_log = []
    while not converged:
        print(f"{timestamp()} BEGINNING POWELL PHASE ITER {counter}")
        h2py.inference.optimize(
            input,
            args.params,
            data["means"],
            data["varcovs"],
            method="powell",
            perturb=args.jitter,
            **kwargs
        )
        print(f"{timestamp()} BEGINNING LBFGSB PHASE ITER {counter}")
        pnames, ps, ll = h2py.inference.optimize(
            input,
            args.params,
            data["means"],
            data["varcovs"],
            method="lbfgsb",
            **kwargs
        )
        delta_ll = ll - ll0
        ll0 = ll
        if args.max_runtime is not None:
            if time.time() - t0 > args.max_runtime:
                converged = True
                print(f"{timestamp()} Optimization exceeded allowed run time")
        if np.abs(delta_ll) < args.threshold:
            converged = True
            print(f"{timestamp()} Optimization converged below criterion")
        counter += 1

