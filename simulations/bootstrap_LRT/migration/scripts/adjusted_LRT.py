
import dpluspy
import pandas
import moments
import numpy as np
import sys


data_file = sys.argv[1]
fit_file0 = sys.argv[2]
param_file0 = sys.argv[3]
fit_file1 = sys.argv[4]
param_file1 = sys.argv[5]
out_file = sys.argv[6]


# Parameters
u = 1.3e-8
weights = (0.5, 0.5)
nested_idx = np.array([4])
# Functions for setting model 0 MLE-based LRT adjustment parameters
func_0 = lambda p0: np.concatenate([p0, [0]])


# Load data
pop_ids, bins, means, varcovs, bootreps = \
    dpluspy.inference.load_bootstrap_reps(data_file)


# Load models
pnames0, params0 = dpluspy.inference._load_params(fit_file0, param_file0)
pnames1, params1 = dpluspy.inference._load_params(fit_file1, param_file1)

model0 = dpluspy.inference.compute_bin_stats(
    fit_file0, sampled_demes=pop_ids, u=u, bins=bins)
model1 = dpluspy.inference.compute_bin_stats(
    fit_file1, sampled_demes=pop_ids, u=u, bins=bins)

ll0 = dpluspy.inference.composite_ll(model0, means, varcovs)
ll1 = dpluspy.inference.composite_ll(model1, means, varcovs)


# Compute the naive LRT statistic
D_naive = 2 * (ll1 - ll0)
p_naive = moments.Godambe.sum_chi2_ppf(D_naive, weights=weights)


# Set up model arguments
_, __, model_args = dpluspy.uncerts.set_up_model_args(
    fit_file1, param_file1, bins=bins, pop_ids=pop_ids, u=u)


# In each case we evaluate with same steps, calculated from model 1 MLE
steps = np.array([params1[i] for i in nested_idx]) * 0.01


# Evaluate `adj0A` at model 0 MLE params, marginalizing nested params
p_lrt0 = func_0(params0)
adj0 = dpluspy.uncerts.LRT_adjust(
    p_lrt0, 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx, 
    steps=steps,
    verbose=0,
)
D_adj0 = adj0 * D_naive
p_adj0 = moments.Godambe.sum_chi2_ppf(D_adj0, weights=weights)


# Evaluate `adj1` at model 1 MLE parameters
adj1 = dpluspy.uncerts.LRT_adjust(
    params1, 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx, 
    steps=steps,
    verbose=0,
)
D_adj1 = adj1 * D_naive
p_adj1 = moments.Godambe.sum_chi2_ppf(D_adj1, weights=weights)


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
pandas.DataFrame(data).to_csv(out_file, index=False)

