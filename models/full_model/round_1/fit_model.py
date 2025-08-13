
import dpluspy
import sys


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
out_file = sys.argv[4]


# Parameters
u = 1.3e-8
perturb = 0.15


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
    method="fmin",
    approx_method="midpoint",
    log=True,
    max_iter=10000,
    verbose=1,
    output=out_file,
    overwrite=True
)
