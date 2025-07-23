
import dpluspy
import numpy as np
import pandas


graph_file =  "den_contraction.yaml"
param_file = "den_contraction_params.yaml"
data_file = "../../../../data/statistics/main/subset_stats.pkl"
out_file = "den_contraction_MLE.yaml"
out_tbl = "contraction_fits.csv"


u = 1.3e-8
pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph=graph_file)


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
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        perturb=0.5,
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
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        method="powell",
        log=True,
        max_iter=100,
        verbose=100,
        output=out_file,
        overwrite=True
    )
    p_log.append(ps)
    if counter > 0:
        diffs = np.abs((np.array(ll_log) - ll) / ll)
        if np.any(diffs < thresh):
            converged = True
        if counter > max_tries:
            raise ValueError("Failed to obtain convergence")
    ll_log.append(ll)
    counter += 1


p_log = np.array(p_log)
trial_log = {"ll": ll_log}
trial_log.update({pname: p_log[:, i] for i, pname in enumerate(pnames)})
pandas.DataFrame(trial_log).to_csv(out_tbl, index=False)


model = dpluspy.inference.compute_bin_stats(out_file, bins=bins, u=u, 
    sampled_demes=pop_ids)
dpluspy.plotting.plot_D_plus_curves(models=model, means=means, varcovs=varcovs,
    bins=bins, pop_ids=pop_ids, out="figure_bestfit.pdf")