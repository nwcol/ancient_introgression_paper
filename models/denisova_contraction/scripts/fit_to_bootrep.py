
import sys
import numpy as np
import pandas
import dpluspy


graph_file = sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
rep_num = int(sys.argv[4])
out_tbl = sys.argv[5]


# PARAMETERS
u = 1.3e-8


pop_ids, bins, _, varcovs, bootreps = dpluspy.inference.load_bootstrap_reps(
    data_file, graph_file)
boot_means = bootreps[rep_num]


ll_log = []
p_log = []
thresh = 0.005
counter = 0
max_tries = 20
converged = False
while not converged:
    dpluspy.inference.optimize(
        graph_file,
        param_file,
        boot_means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        perturb=0.1,
        method="lbfgsb",
        log=True,
        max_iter=2000,
        verbose=100,
        output="temp.yaml",
        overwrite=True
    )
    pnames, ps, ll = dpluspy.inference.optimize(
        "temp.yaml",
        param_file,
        boot_means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        method="powell",
        log=True,
        max_iter=100,
        verbose=100,
    )
    if counter > 0:
        diffs = np.abs((np.array(ll_log) - ll) / ll)
        if np.any(diffs < thresh):
            converged = True
        if counter > max_tries:
            raise ValueError("Failed to obtain convergence")
    ll_log.append(ll)
    counter += 1


p_log = np.array(p_log)
trial_log = {"ll": [ll]}
trial_log.update({pname: [ps[i]] for i, pname in enumerate(pnames)})
pandas.DataFrame(trial_log).to_csv(out_tbl, index=False)
