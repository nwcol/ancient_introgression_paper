
import demes
import numpy as np
import moments
from h2py.utils import timestamp
from moments.Demes.Inference import (
    _update_builder,
    _get_demes_dict,
    _get_params_dict,
)
import pickle
import argparse
import pandas
import time


def get_args():
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


def print_params(pnames, ps, ll):
    print(timestamp(), f"ll={ll:.6}")
    for pname, p in zip(pnames, ps):
        print(f"{pname}\t{p:.4}")
    return


if __name__ == "__main__":
    args = get_args()

    with open(args.data, "rb") as fin:
        data = pickle.load(fin)

    L = data["L"]
    fs = data["fs"]
    u = args.u
    uL = u * L

    t0 = time.time()
    converged = False

    kwargs = {
        "uL": uL,
        "verbose": args.report,
        "output": args.out,
        "overwrite": True,
        "maxiter": args.max_iter
    }

    print(f"{timestamp()} BEGINNING INITIAL FMIN PHASE")
    pnames, ps, ll0 = moments.Demes.Inference.optimize(
        args.graph,
        args.params,
        fs,
        perturb=args.perturb,
        method="fmin",
        **kwargs
    )
    print_params(pnames, ps, ll0)

    input = args.out
    converged = False
    counter = 0
    while not converged:
        print(f"{timestamp()} BEGINNING POWELL PHASE ITER {counter}")
        pnames, ps, ll = moments.Demes.Inference.optimize(
            input,
            args.params,
            fs,
            method="powell",
            perturb=args.jitter,
            **kwargs
        )
        print_params(pnames, ps, ll)

        print(f"{timestamp()} BEGINNING LBFGSB PHASE ITER {counter}")
        pnames, ps, ll = moments.Demes.Inference.optimize(
            input,
            args.params,
            fs,
            method="lbfgsb",
            **kwargs
        )
        print_params(pnames, ps, ll)

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

