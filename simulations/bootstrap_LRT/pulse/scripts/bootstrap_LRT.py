# Fit a null and complex model to bootstrap replicates

import dpluspy
import pandas
import numpy as np
import sys


data_file = sys.argv[1]
graph_file0 = sys.argv[2]
param_file0 = sys.argv[3]
graph_file1 = sys.argv[4]
param_file1 = sys.argv[5]
out_file = sys.argv[6]


# Parameters
u = 1.3e-8


pop_ids, bins, means, varcovs, bootreps = \
    dpluspy.inference.load_bootstrap_reps(data_file)


def fit_model(graph_file, param_file, bootrep):
    """
    Fit a model, forcing convergence to within 1% LL.
    """
    ll_log = []
    thresh = 0.01
    counter = 0
    max_tries = 5
    converged = False
    while not converged:
        if counter > max_tries:
            break
        ret = dpluspy.inference.optimize(
            graph_file,
            param_file,
            bootrep,
            varcovs,
            pop_ids=pop_ids,
            bins=bins,
            u=u,
            perturb=0.05,
            method="powell",
            log=True,
            max_iter=100,
            verbose=200
        )
        ll = ret[-1]
        if counter > 0:
            diffs = np.abs((np.array(ll_log) - ll) / ll)
            if np.any(diffs < thresh):
                converged = True
        ll_log.append(ll)
        counter += 1
    return ret
 

rep_ll0s = []
rep_ll1s = []
for ii, bootrep in enumerate(bootreps):
    _, p0, ll0 = fit_model(graph_file0, param_file0, bootrep)
    _, p1, ll1 = fit_model(graph_file1, param_file1, bootrep)
    rep_ll0s.append(ll0)
    rep_ll1s.append(ll1)

data = {
    "ll0": rep_ll0s,
    "ll1": rep_ll1s
}

pandas.DataFrame(data).to_csv(out_file, index=False)


