
import demes
import moments
import numpy as np
import dpluspy
import sys
import pandas


# Arguments
model_graph = sys.argv[1]
model_params = sys.argv[2]
data_fname = sys.argv[3]
out_fstem = sys.argv[4]

log_fname = f"{out_fstem}_log.csv"
out_fname = f"{out_fstem}.csv"


# Constant parameters
u = 1.3e-08
num_tries = 5
max_cycles = 3
max_iter_lbfgsb = 100
max_iter_powell = 100
perturb = 0.10
jitter = 0.03
threshold = 1e-3


pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_fname, graph=model_graph)

# Keyword arguments for optimization
kwargs = {
    "pop_ids": pop_ids,
    "bins": bins,
    "u": u,
    "log": True,
    "verbose": 10000
}

log = []


def make_df(pnames, params, ll):
    mapping = {x: [y] for x, y in zip(["ll"] + pnames, [ll] + list(params))}
    return pandas.DataFrame(mapping)


def fit_model(graph_fname, param_fname, log_id="", init_perturb=0):
    last_ll = None
    counter = 0
    converged = False
    while not converged:
        print(f"INITIATING CYCLE {counter}")
        if counter == 0:
            _init_perturb = init_perturb
        else:
            _init_perturb = jitter
        dpluspy.inference.optimize(
            graph_fname,
            param_fname,
            means,
            varcovs,
            method="lbfgsb",
            max_iter=max_iter_lbfgsb,
            perturb=_init_perturb,
            output="temp.yaml",
            overwrite=True,
            **kwargs
        )
        pnames, params, ll = dpluspy.inference.optimize(
            "temp.yaml",
            param_fname,
            means,
            varcovs,
            method="powell",
            max_iter=max_iter_powell,
            perturb=0,
            output="temp.yaml",
            overwrite=True,
            **kwargs
        )
        if counter > 0:
            diff = np.abs(last_ll - ll)
            if diff < threshold:
                converged = True
            if counter > max_cycles:
                converged = True
        df = make_df(
            ["trial"] + pnames, [f"{log_id}_{counter}"] + list(params), ll)
        log.append(df)
        last_ll = ll
        counter += 1
        print(f"{log_id}_{counter} fitted")
    return pnames, params, ll


fitted_models = {}

for ii in range(num_tries):
    print(f"INITIATING TRY {ii}")
    pnames, params, ll = fit_model(
        model_graph, 
        model_params,
        f"try_{ii}",
        init_perturb=perturb
    )
    fitted_models[ii] = (ll, params)
    lls_so_far = [fitted_models[x][0] for x in fitted_models]
    if np.any(np.abs(lls_so_far[:-1] - ll) < threshold):
        break

lls = [fitted_models[x][0] for x in fitted_models]
mle_idx = np.argmax(lls)
max_ll, mle_params = fitted_models[mle_idx]
results = make_df(pnames, list(mle_params), max_ll)
results.to_csv(out_fname, index=False)
pandas.concat(log).to_csv(log_fname, index=False)
