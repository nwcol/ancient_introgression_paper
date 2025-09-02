
import dpluspy
import sys


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
out_file = sys.argv[4]


# Parameters
u = 1.3e-8

# Absolute difference in ll units required for cessation of model fitting
threshold = 1e-5
max_tries = 10
perturb = 0.7
verbose = 1000


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
    log=True,
    max_iter=10000,
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
