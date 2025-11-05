
import demes
import moments
import numpy as np
import dpluspy
import sys
import pandas


# Arguments
model_0_graph = sys.argv[1]
model_0_params = sys.argv[2]
model_1_graph = sys.argv[3]
model_1_params = sys.argv[4]
data_fname = sys.argv[5]
out_fstem = sys.argv[6]

log_fname = f"{out_fstem}_log.csv"
out_fname = f"{out_fstem}.csv"


# Key parameter check
if "sd" in model_0_graph:
    pulse_pname = "p_S_D"
elif "hn" in model_0_graph:
    pulse_pname = "p_H_N"
else:
    pulse_pname = ""


# Constant parameters
u = 1.3e-08
num_tries = 5
max_cycles = 5
perturb = 0.30
jitter = 0.05
threshold = 1e-3


pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_fname, graph=model_0_graph)

# Keyword arguments for optimization
kwargs = {
    "pop_ids": pop_ids,
    "bins": bins,
    "u": u,
    "log": True,
    "verbose": 10
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
        if counter == 0:
            perturb_fmin = init_perturb
        else:
            perturb_fmin = jitter
        dpluspy.inference.optimize(
            graph_fname,
            param_fname,
            means,
            varcovs,
            method="fmin",
            max_iter=1000,
            perturb=perturb_fmin,
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
            max_iter=100,
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


def transpose_params(
    pnames0, 
    params0, 
    graph_fname1, 
    param_fname1,
    output_fname,
    mapping={}
):
    """
    Set the graph of an alternate model to the MLE parameters of a null model
    using the null model parameterization, and save the result.
    """
    builder = moments.Demes.Inference._get_demes_dict(graph_fname1)
    if param_fname1 is not None:
        options = moments.Demes.Inference._get_params_dict(param_fname1)
        pnames1, params1 = dpluspy.inference._load_params(
            graph_fname1, param_fname1)
        for i, pname in enumerate(pnames1): 
            if pname in mapping:
                params1[i] = mapping[pname]
            elif pname in pnames0:
                idx = pnames0.index(pname)
                params1[i] = params0[idx]
        params = params1
    builder = moments.Demes.Inference._update_builder(builder, options, params)
    graph = demes.Graph.fromdict(builder)
    demes.dump(graph, output_fname)
    return


# fit null models
null_models = {}

for ii in range(num_tries):
    pnames0, params0, ll0 = fit_model(
        model_0_graph, 
        model_0_params,
        f"null_{ii}",
        init_perturb=perturb
    )
    null_models[ii] = (ll0, params0)
    lls_so_far = [null_models[x][0] for x in null_models]
    if np.any(np.abs(lls_so_far[:-1] - ll0) < threshold):
        break

lls0 = [null_models[x][0] for x in null_models]
mle_idx0 = np.argmax(lls0)
max_ll0, mle_params0 = null_models[mle_idx0]


alt_models = {}

# fit an alternative model from ~model 0 MLE
transpose_params(
    pnames0, 
    mle_params0, 
    model_1_graph, 
    model_1_params, 
    "init1.yaml", 
    mapping={pulse_pname: 1e-5}
)
pnames1, params1, ll1 = fit_model(
    "init1.yaml", model_1_params, "alt_0", init_perturb=0)
alt_models[0] = (ll1, params1)


transpose_params(
    pnames0, 
    mle_params0, 
    model_1_graph, 
    model_1_params, 
    "init2.yaml", 
    mapping={pulse_pname: 1e-2}
)
for ii in range(1, num_tries):
    pnames1, params1, ll1 = fit_model(
        "init2.yaml", 
        model_1_params,
        f"alt_{ii}",
        init_perturb=perturb
    )
    alt_models[ii] = (ll1, params1)
    lls_so_far = [alt_models[x][0] for x in alt_models]
    diffs = np.abs(lls_so_far[:-1] - ll1)
    if np.any(diffs < threshold) and ll1 > max_ll0:
        break

lls1 = [alt_models[x][0] for x in alt_models]
mle_idx1 = np.argmax(lls1)
max_ll1, mle_params1 = alt_models[mle_idx1]

results = pandas.concat(
    [make_df(["model"] + pnames0, ["null"] + list(mle_params0), max_ll0),
    make_df(["model"] + pnames1, ["alt"] + list(mle_params1), max_ll1)])

results.to_csv(out_fname, index=False)
pandas.concat(log).to_csv(log_fname, index=False)
