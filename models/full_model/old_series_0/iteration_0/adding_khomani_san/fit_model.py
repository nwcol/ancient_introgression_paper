
import dpluspy
import numpy as np
import pandas


graph_file =  "basal_eurasian_model_vin_Bherer_7.yaml"
param_file = "addition_params.yaml"
data_file = "../../../../data/statistics/main/subset_stats.pkl"
out_file = "addition_fit.yaml"

u = 1.3e-8
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
    method="lbfgsb",
    log=True,
    max_iter=500,
    verbose=10,
    output="temp.yaml",
    overwrite=True
)
pnames, ps, ll = dpluspy.inference.optimize(
    "temp.yaml",
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    method="powell",
    log=True,
    max_iter=100,
    verbose=10,
    output=out_file,
    overwrite=True
)

model = dpluspy.inference.compute_bin_stats(out_file, bins=bins, u=u, 
    sampled_demes=pop_ids)
dpluspy.plotting.plot_D_plus_curves(models=model, means=means, varcovs=varcovs,
    bins=bins, pop_ids=pop_ids, out="best_fit.pdf")