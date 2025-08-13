
import dpluspy
import numpy as np
import pandas


graph_file =  "model.yaml"
param_file = "model_params.yaml"
data_file = "../../../../data/statistics/main/subset_stats.pkl"
out_file = "model_fit.yaml"


u = 1.3e-8
pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph=graph_file)

pnames, ps, ll = dpluspy.inference.optimize(
    graph_file,
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    method="lbfgsb",
    log=True,
    max_iter=4000,
    verbose=10,
    output=out_file,
    overwrite=True
)


model = dpluspy.inference.compute_bin_stats(out_file, bins=bins, u=u, 
    sampled_demes=pop_ids)
dpluspy.plotting.plot_D_plus_curves(models=model, means=means, varcovs=varcovs,
    bins=bins, pop_ids=pop_ids, out="fit_figure.pdf")