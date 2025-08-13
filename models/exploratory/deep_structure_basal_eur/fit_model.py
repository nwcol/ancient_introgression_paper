
import demes
import dpluspy
import numpy as np
import pandas
import moments
import sys


graph_file = "model.yaml"
param_file = "model_params.yaml"
data_file = "../../../analyses/recombination_maps/stats/subset_Bherer.pkl"
out_file = "fitted_model.yaml"


# Parameters
u = 1.3e-8

# Absolute difference in ll units required for cessation of model fitting
threshold = 1e-4
max_tries = 10
perturb = 0
verbose = 20


pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph=graph_file)


dpluspy.inference.optimize(
    graph_file,
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    perturb=perturb,
    method="lbfgsb",
    log=True,
    max_iter=1000,
    verbose=verbose,
    output=out_file,
    overwrite=True
)


# Fit the model to convergence using the Powell algorithm
converged = False
counter = 0
ll_log = []
while not converged:
    pnames, ps, ll = dpluspy.inference.optimize(
        out_file,
        param_file,
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        method="powell",
        log=True,
        max_iter=100,
        verbose=verbose,
        output=out_file,
        overwrite=True
    )
    if counter > 0:
        last_ll = ll_log[-1]
        ll_diff = ll - last_ll
        if ll_diff < threshold:
            converged = True
        if counter > max_tries:
            converged = True
    ll_log.append(ll)
    counter += 1


model = dpluspy.inference.compute_bin_stats(
    out_file, bins=bins, u=u, sampled_demes=pop_ids)
dpluspy.plotting.plot_D_plus_curves(models=model, means=means, varcovs=varcovs,
    bins=bins, pop_ids=pop_ids, out="fig_model_fit.pdf")