# Fit a model to data

import dpluspy
import numpy as np
import sys


data_file = sys.argv[1]
graph_file = sys.argv[2]
param_file = sys.argv[3]
out_file = sys.argv[4]


# Parameters
u = 1.3e-8


pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(data_file)


def fit_model():
    """
    Fit a model, forcing convergence to within 1% LL.
    """
    ll_log = []
    thresh = 0.01
    counter = 0
    max_tries = 20
    converged = False
    while not converged:
        ret = dpluspy.inference.optimize(
            graph_file,
            param_file,
            means,
            varcovs,
            pop_ids=pop_ids,
            bins=bins,
            u=u,
            perturb=0.05,
            method="powell",
            log=True,
            max_iter=100,
            verbose=200,
            output=out_file,
            overwrite=True
        )
        ll = ret[-1]
        if counter > 0:
            diffs = np.abs((np.array(ll_log) - ll) / ll)
            if np.any(diffs < thresh):
                converged = True
            if counter > max_tries:
                raise ValueError("Failed to obtain convergence")
        ll_log.append(ll)
        counter += 1
    return ret
 

fit_model()

