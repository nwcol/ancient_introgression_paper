"""
Fit nested null/alternative models to the SFS. Alternative model fits begin
at the MLE of the null model. Parameters that do not exist in the null model
start at values specified by the  --nested argument.
"""

import argparse
import demes
from h2py.utils import timestamp
import moments
from moments.Demes.Inference import (
    _update_builder,
    _get_demes_dict,
    _get_params_dict,
    _set_up_params_and_bounds,
)
import numpy as np
import pickle


def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--data", required=True)
    parser.add_argument("-g0", "--graph0", required=True)
    parser.add_argument("-p0", "--params0", required=True)
    parser.add_argument("-g1", "--graph1", required=True)
    parser.add_argument("-p1", "--params1", required=True)
    parser.add_argument("--nested", nargs="*",
                        help="nested parameter-value pairs, e.g. T_x 1000 p_x 0.0")
    parser.add_argument("-o", "--out", required=True,
                        help="prefix for output files, e.g., {out}_0.yaml")
    # Optional arguments
    parser.add_argument("-u", type=float, default=1.3e-8)
    parser.add_argument("--max-iter", type=int, default=1000)
    parser.add_argument("--perturb", type=float, default=0)
    return parser.parse_args()


if __name__ == "__main__":
    args = get_args()
    # Load input data
    with open(args.data, "rb") as fin:
        data = pickle.load(fin)
    # Set up input data
    uL = float(args.u * data["L"])
    sfs = data["sfs"]

    # Set up a function for fitting models
    def fit_model(graph, params, perturb, out):
        # Invariant arguments
        kwargs = {
            "maxiter": args.max_iter,
            "output": out,
            "overwrite": True,
            "uL": uL,
            "verbose": 100,
        }
        # Stage 1: LBFGSB
        pnames, ps, ll = moments.Demes.Inference.optimize(
            graph,
            params,
            sfs,
            method="lbfgsb",
            perturb=perturb,
            **kwargs
        )
        moments.Demes.Inference._counter = 0
        print(timestamp(), f"{graph} lbfgsb stage finished, ll={ll}")
        # Stage 2: powell
        pnames, ps, ll = moments.Demes.Inference.optimize(
            graph,
            params,
            sfs,
            method="powell",
            **kwargs
        )
        moments.Demes.Inference._counter = 0
        print(timestamp(), f"{graph} fmin stage finished, ll={ll}")
        # Add optimization information to the graph's metadata
        g = demes.load(out)
        g.metadata["opt_info"] = {"ll": -ll, "u": args.u}
        demes.dump(g, out)
        return pnames, ps, -ll

    # Fit the null model to convergence
    pnames0, ps0, ll0 = fit_model(args.graph0, args.params0, args.perturb,
                                  f"{args.out}_H0.yaml")
    print(timestamp(), f"finished fitting model {args.graph0}!")
    # Get default values for parameters absent from the null model
    nest_params = {key: float(val)
                   for key, val in zip(args.nested[::2], args.nested[1::2])}
    print(timestamp(), f"nesting: {nest_params}")
    builder1 = _get_demes_dict(args.graph1)
    options1 = _get_params_dict(args.params1)
    pnames1, ps1 = _set_up_params_and_bounds(options1, builder1)[:2]
    # Set up the alternative model at MLE parameters of the null model
    ps1_work = []
    for pname1 in pnames1:
        if pname1 in pnames0:
            idx = pnames0.index(pname1)
            ps1_work.append(ps0[idx])
        else:
            assert pname1 in nest_params
            val = nest_params[pname1]
            ps1_work.append(val)
    builder1 = _update_builder(builder1, options1, ps1_work)
    work_graph1 = demes.Graph.fromdict(builder1)
    # Write the constructed alternative model to disk
    model1_start = f"{args.out}_H1_init.yaml"
    demes.dump(work_graph1, model1_start)
    # Fit the null model to convergence
    perturb1 = 0
    pnames1, pvals1, ll1 = fit_model(model1_start, args.params1, perturb1,
                                     f"{args.out}_H1.yaml")
    print(timestamp(), f"finished fitting model {args.graph1}!")
