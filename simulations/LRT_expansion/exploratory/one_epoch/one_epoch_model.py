
import dpluspy
import demes 
import moments
import msprime
import numpy as np
import pandas
import pickle
import sys


rep = sys.argv[1]
out_file = f"one_epoch_results_{rep}.csv"


# Define models and parameters
model0_file = "models/graph0.yaml"
params0_file = "models/params0.yaml"
model1_file = "models/graph1.yaml"
params1_file = "models/params1.yaml"


# Parameters
r = 1e-8
u = 1.3e-8
u_std = 1e-9
L = int(5e6)
n_reps = 100
pop_ids = ["popX"]
samples = {pop_id: 1 for pop_id in pop_ids}
bins = np.logspace(-6, -2, 17)


def simulate():
    graph = demes.load(model0_file)
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
    intervals = [[1, L + 1, L + 1]]
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

replicates = dpluspy.bootstrapping.get_bootstrap_reps(
    stats, num_reps=100, weighted=False)
varcovs = dpluspy.bootstrapping.compute_varcovs(replicates)
means = dpluspy.bootstrapping.means_across_regions(stats)
boot_stats = {"means": means, "varcovs": varcovs, "replicates": replicates,
    "bins": bins, "pop_ids": pop_ids}

# Save result
with open("statistics.pkl", "wb") as fout:
    pickle.dump(boot_stats, fout)


def fit_model(graph_file, param_file):
    """
    Fit a model in two stages; return fitted parameters without writing output
    to disk
    """
    dpluspy.inference.optimize(
        graph_file,
        param_file,
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        perturb=0.2,
        method="lbfgsb",
        log=True,
        max_iter=10000,
        verbose=100,
        output="intermediate.yaml",
        overwrite=True
    )
    ret = dpluspy.inference.optimize(
        "intermediate.yaml",
        param_file,
        means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        method="powell",
        log=True,
        max_iter=100,
        verbose=100
    )
    return ret
 

# Fit the underlying and complex models; try to force ll1 > ll0
pnames0, params0, ll0 = fit_model(model0_file, params0_file)
ll1 = -1e10
max_tries = 20
tries = 0
while ll1 <= ll0:
    pnames1, params1, ll1 = fit_model(model1_file, params1_file)
    tries += 1
    if tries > max_tries:
        raise ValueError("Cannot achieve ll1 > ll0")

# We have one nested parameter, which is not at a boundary
chi_sqr_weights = (0.5, 0.5)
D_naive = 2 * (ll1 - ll0)
p_naive = moments.Godambe.sum_chi2_ppf(D_naive, weights=chi_sqr_weights)

# Set up model arguments
_, __, model_args = dpluspy.uncerts.set_up_model_args(
    model1_file, params1_file, bins=bins, pop_ids=pop_ids, u=u)

# Compute LRT adjustment at the simple model MLE parameter values:
p0 = np.array([params0[0], 1, params0[0]])
steps = np.array([params1[1], params0[0]]) * 0.01
nested_idx = np.array([1, 2])

adj_simple = dpluspy.uncerts.LRT_adjust(
    p0, 
    model_args, 
    means, 
    varcovs, 
    replicates, 
    nested_idx, 
    verbose=True,
    steps=steps
)

# Compute the LR test statistic
D_simple = adj_simple * D_naive
p_simple = moments.Godambe.sum_chi2_ppf(D_simple, weights=chi_sqr_weights)

# Now compute LRT adjustment at the complex model MLE parameter values,
p1 = np.array([params1[0], params1[1], params1[0]])
steps = np.array([params1[1], params1[0]]) * 0.01

adj_complex = dpluspy.uncerts.LRT_adjust(
    p1, 
    model_args, 
    means, 
    varcovs, 
    replicates, 
    nested_idx, 
    verbose=True
)

D_complex = adj_complex * D_naive
p_complex = moments.Godambe.sum_chi2_ppf(D_complex, weights=chi_sqr_weights)

data = {
    "ll0": [ll0],
    "ll1": [ll1],
    "adj_simple": [adj_simple],
    "adj_complex": [adj_complex],
    "D_naive": [D_naive],
    "p_naive": [p_naive],
    "D_simple": [D_simple],
    "p_simple": [p_simple],
    "D_complex": [D_complex],
    "p_complex": [p_complex]
}
data.update({f"{pnames0[i]}_simple": [params0[i]] for i in range(len(params0))})
data.update({f"{pnames1[i]}_complex": [params1[i]] for i in range(len(params1))})
pandas.DataFrame(data).to_csv(out_file, index=False)
