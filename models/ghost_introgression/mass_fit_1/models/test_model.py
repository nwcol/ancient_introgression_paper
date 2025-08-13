
import demes
import dpluspy
import numpy as np
import pandas
import moments
import sys


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = "../../../analyses/recombination_maps/stats/subset_Bherer.pkl"
out_file = f"scratch_{graph_file}"


# Parameters
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
    perturb=0,
    method="fmin",
    approx_method="midpoint",
    log=True,
    max_iter=35,
    verbose=10,
    output=out_file,
    overwrite=True
)

