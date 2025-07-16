
import dpluspy
import demes 
import moments
import msprime
import numpy as np
import pandas
import pickle
import sys


rep = sys.argv[1]
model_idx = sys.argv[2]


# To be explicit about the underlying model used in simulation:
true_model = f"models/graph{model_idx}.yaml"


out_file = f"symmetry_LRT_model{model_idx}_rep{rep}.csv"
stats_file = f"symmetry_LRT_model{model_idx}_rep{rep}_stats.pkl"
out_graph0 = f"symmetry_LRT_model{model_idx}_rep{rep}_fit0.yaml"
out_graph1 = f"symmetry_LRT_model{model_idx}_rep{rep}_fit1.yaml"


# Define models and parameters
model0_file = "models/graph0.yaml"
params0_file = "models/params0.yaml"
model1_file = "models/graph1.yaml"
params1_file = "models/params1.yaml"


# Parameters
r = 1e-8
u = 1.3e-8
u_std = 2e-9
L = int(5e6)
n_reps = 100
pop_ids = ["popX", "popY"]
samples = {pop_id: 1 for pop_id in pop_ids}
bins = np.logspace(-6, -2, 17)


def simulate():
    graph = demes.load(true_model)
    demog = msprime.Demography.from_demes(graph)
    tss = msprime.sim_ancestry(
        demography=demog,
        samples=samples,
        sequence_length=L,
        recombination_rate=r,
        num_replicates=n_reps,
    )
    mtss = [msprime.sim_mutations(ts, rate=np.random.normal(u, u_std)) 
        for ts in tss]
    return mtss


def parse_stats(mtss):
    # Parse statistics from simulation results
    intervals = [[1, L + 1]]
    pop_mapping = {x: [x] for x in samples} 
    stats = dict()
    for ii, ts in enumerate(mtss):
        if ii == 0:
            stats_ii = dpluspy.parsing.parse_stats(
                ts, 
                get_denoms=True,
                pop_mapping=pop_mapping, 
                r=r, 
                r_bins=bins, 
                intervals=intervals,
                chrom=ii, 
                ts_sample_ids=pop_ids, 
                overhang=None
            )
        else:
            stats_ii = dpluspy.parsing.parse_stats(
                ts, 
                get_denoms=False,
                pop_mapping=pop_mapping, 
                r=r, 
                r_bins=bins, 
                intervals=intervals,
                chrom=ii, 
                ts_sample_ids=pop_ids, 
                overhang=None
            )
            for jj in range(len(intervals)):
                stats_ii[(ii, jj)]["denoms"] = stats[(0, jj)]["denoms"]
        stats.update(stats_ii)
    return stats


mtss = simulate()
stats = parse_stats(mtss)

bootreps = dpluspy.bootstrapping.get_bootstrap_reps(
    stats, num_reps=100, weighted=False)
varcovs = dpluspy.bootstrapping.compute_varcovs(bootreps)
means = dpluspy.bootstrapping.means_across_regions(stats)

# Save result
boot_stats = {
    "means": means, 
    "varcovs": varcovs, 
    "replicates": bootreps,
    "bins": bins, 
    "pop_ids": pop_ids
}
with open(stats_file, "wb") as fout:
    pickle.dump(boot_stats, fout)


def fit_model(graph_file, param_file, output):
    """
    Fit a model in two stages, forcing convergence to within 0.1% LL.
    """
    ll_log = []
    thresh = 0.001
    counter = 0
    max_tries = 20
    converged = False
    while not converged:
        _ = dpluspy.inference.optimize(
            graph_file,
            param_file,
            means,
            varcovs,
            pop_ids=pop_ids,
            bins=bins,
            u=u,
            perturb=0.3,
            method="lbfgsb",
            log=True,
            verbose=0,
            max_iter=2000,
            output="round1.yaml",
            overwrite=True
        )
        ret = dpluspy.inference.optimize(
            "round1.yaml",
            param_file,
            means,
            varcovs,
            pop_ids=pop_ids,
            bins=bins,
            u=u,
            method="powell",
            log=True,
            max_iter=100,
            output=output,
            verbose=1000,
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
 

# Fit the underlying and complex models; try to force ll1 > ll0
pnames0, params0, ll0 = fit_model(model0_file, params0_file, out_graph0)
ll1 = -1e10
max_tries = 20
tries = 0
while ll1 <= ll0:
    pnames1, params1, ll1 = fit_model(model1_file, params1_file, out_graph1)
    tries += 1
    if tries > max_tries:
        raise ValueError("Cannot achieve ll1 > ll0")


# We have one nested parameter which is not at the boundary
chi_sqr_weights = (0, 1)
D_naive = 2 * (ll1 - ll0)
p_naive = moments.Godambe.sum_chi2_ppf(D_naive, weights=chi_sqr_weights)

# Set up model arguments
_, __, model_args = dpluspy.uncerts.set_up_model_args(
    model1_file, 
    params1_file, 
    bins=bins, 
    pop_ids=pop_ids,
    u=u
)
nested_idx = np.array([3])


# Evaluate at model 0 MLE parameters, setting N_Y = N_X_Y
p_lrt0 = np.concatenate([params0, [params0[-1]]])
steps0 = np.array([p_lrt0[-1] * 0.01])
adj0 = dpluspy.uncerts.LRT_adjust(
    p_lrt0, 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx, 
    steps=steps0,
    verbose=0,
)
D_adj0 = adj0 * D_naive
p_adj0 = moments.Godambe.sum_chi2_ppf(D_adj0, weights=chi_sqr_weights)


# Evaluate at model 1 MLE parameters
p_lrt1 = params1
steps1 = np.array([p_lrt1[-1] * 0.01])
adj1 = dpluspy.uncerts.LRT_adjust(
    p_lrt1, 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx, 
    steps=steps1,
    verbose=0,
)
D_adj1 = adj1 * D_naive
p_adj1 = moments.Godambe.sum_chi2_ppf(D_adj1, weights=chi_sqr_weights)


data = {
    "ll0": [ll0],
    "ll1": [ll1],
    "D_naive": [D_naive],
    "p_naive": [p_naive],
    "adj0": [adj0],
    "D_adj0": [D_adj0],
    "p_adj0": [p_adj0],
    "adj1": [adj1],
    "D_adj1": [D_adj1],
    "p_adj1": [p_adj1],
}
data.update({f"{pnames0[i]}_simple": [params0[i]] for i in range(len(params0))})
data.update({f"{pnames1[i]}_complex": [params1[i]] for i in range(len(params1))})
pandas.DataFrame(data).to_csv(out_file, index=False)
