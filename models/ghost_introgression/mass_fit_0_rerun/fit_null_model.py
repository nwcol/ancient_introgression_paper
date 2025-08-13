
import demes
import dpluspy
import numpy as np
import pandas
import moments
import sys


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
out_file = sys.argv[4]
out_tbl = sys.argv[5]


# Parameters
thresh = 0.01
counter = 0
max_tries = 10
converged = False
perturb = 0.50


u = 1.3e-8
pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph=graph_file)


ll_log = []
p_log = []
while not converged:
    dpluspy.inference.optimize(
        graph_file,
        param_file,
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        perturb=perturb,
        method="lbfgsb",
        log=True,
        max_iter=1000,
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
    p_log.append(ps)
    if counter > 0:
        diff = np.abs((np.max(ll_log) - ll) / ll)
        if diff < thresh:
            converged = True
        if counter > max_tries:
            converged = True
    ll_log.append(ll)
    counter += 1


# Write output YAML file
opt_idx = int(np.argmax(ll_log))
opt_params = p_log[opt_idx]
metadata = {"opt_info": {"ll": ll_log[opt_idx]}}
builder = moments.Demes.Inference._get_demes_dict(graph_file)
options = moments.Demes.Inference._get_params_dict(param_file)    
builder = moments.Demes.Inference._update_builder(builder, options, opt_params)
g = demes.Graph.fromdict(builder)
g.metadata = metadata
demes.dump(g, out_file)


p_log = np.array(p_log)
trial_log = {"ll": ll_log}
trial_log.update({pname: p_log[:, i] for i, pname in enumerate(pnames)})
pandas.DataFrame(trial_log).to_csv(out_tbl, index=False)
