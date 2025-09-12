
import dpluspy
import sys
import time
import numpy as np


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
out_file = sys.argv[4]


# Parameters
u = 1.3e-8
perturb = 1
verbose = 10
criterion = 1e-4
max_run_time = 10.0


pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph=graph_file)


t0 = time.time()
converged = False


_, __, ll0 = dpluspy.inference.optimize(
    graph_file,
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    perturb=perturb,
    method="lbfgsb",
    approx_method="midpoint",
    log=True,
    max_iter=1,
    verbose=verbose,
    output=out_file,
    overwrite=True
)

while not converged:
    _, __, ll = dpluspy.inference.optimize(
        out_file,
        param_file,
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        method="fmin",
        approx_method="midpoint",
        log=True,
        max_iter=50,
        verbose=verbose,
        output=out_file,
        overwrite=True
    )
    delta_ll = ll - ll0
    ll0 = ll
    if time.time() - t0 > max_run_time:
        converged = True 
        print("Optimization exceeded allowed run time")
    if np.abs(delta_ll) < criterion:
        converged = True
        print("Optimization converged below criterion")
