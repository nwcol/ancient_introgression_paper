"""
Model 0 adjustment strategies: 
    0_0: marginal model 0 parameters. T_0 ~ 0, N_1 = N_0
    0_1: T_0 fixed at model 1 MLE, N_1 = N_0
    0_2: same as 0_1 but treats T_0 as though fixed or not nested
"""

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
weights = (0, 0.5, 0.5)
nested_idx = np.array([1, 2])
nested_idx_2 = np.array([2])
# Functions for setting model 0 MLE-based LRT adjustment parameters
param_func0 = lambda p0: np.concatenate([p0, [1e-10, p0[0]]])
step_func0 = lambda p0: params1[1:] * 0.01
param_func1 = lambda p0: np.concatenate([p0, [params1[1], p0[0]]])
step_func2 = lambda p0: params1[2:] * 0.01


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


# Evaluate LRT adjustments
try:
    adj0_0 = dpluspy.uncerts.LRT_adjust(
        param_func0(params0), 
        model_args, 
        means, 
        varcovs, 
        bootreps, 
        nested_idx, 
        steps=step_func0(params0),
        verbose=0,
    )
except:
    adj0_0 = np.nan


adj0_1 = dpluspy.uncerts.LRT_adjust(
    param_func1(params0), 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx, 
    steps=step_func0(params0),
    verbose=0,
)


adj0_2 = dpluspy.uncerts.LRT_adjust(
    param_func1(params0), 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx_2, 
    steps=step_func2(params0),
    verbose=0,
)


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


data = {
    "ll0": [ll0],
    "ll1": [ll1],
    "adj0_0": [adj0_0],
    "adj0_1": [adj0_1],
    "adj0_2": [adj0_2],
    "adj1": [adj1],
}
pandas.DataFrame(data).to_csv(out_file, index=False)

