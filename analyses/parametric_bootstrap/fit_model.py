
import h2py
import sys
import time
import numpy as np
import pickle


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
out_file = sys.argv[4]


# Parameters
u = 1.3e-8
perturb = 0.05
jitter = 0.02
verbose = 200
threshold = 1e-2
max_run_time = 86400
max_iter_primary = 1000
max_iter_secondary = 100

with open(data_file, "rb") as fin:
    data = pickle.load(fin)

t0 = time.time()
converged = False

kwargs = {
    "pops": data["pops"],
    "r_bins": data["bins"],
    "u": u,
    "use_log": True,
    "report": verbose,
    "output": out_file,
    "overwrite": True
}


print("INITIAL FMIN PHASE")
_, __, ll0 = h2py.inference.optimize(
    graph_file,
    param_file,
    data["means"],
    data["varcovs"],
    method="fmin",
    max_iter=max_iter_primary,
    perturb=perturb,
    **kwargs
)


# Fit the model to convergence, alternating between Powell and LBFGSB algorithms
converged = False
counter = 0
ll_log = []
while not converged:
    print(f"POWELL PHASE ITER {counter}")
    h2py.inference.optimize(
        out_file,
        param_file,
        data["means"],
        data["varcovs"],
        method="powell",
        max_iter=max_iter_secondary,
        perturb=jitter,
        **kwargs
    )
    print(f"LBFGSB PHASE ITER {counter}")
    pnames, ps, ll = h2py.inference.optimize(
        out_file,
        param_file,
        data["means"],
        data["varcovs"],
        method="lbfgsb",
        max_iter=max_iter_secondary,
        **kwargs
    )
    delta_ll = ll - ll0
    ll0 = ll
    if time.time() - t0 > max_run_time:
        converged = True 
        print("Optimization exceeded allowed run time")
    if np.abs(delta_ll) < threshold:
        converged = True
        print("Optimization converged below criterion")
    counter += 1

