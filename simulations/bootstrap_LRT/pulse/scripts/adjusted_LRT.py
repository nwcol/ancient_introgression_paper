"""
Model 0 adjustment strategies: 
    0_0: marginal model 0 parameters. p_XY = 0, T_pXY = T_XY
    0_1: T_pXY fixed at model 1 MLE, p_XY = 0
    0_2: same as 0_1 but treats T_pXY as though fixed or not nested
"""

import dpluspy
import pandas
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
nested_idx = np.array([4, 5])
nested_idx_2 = np.array([5])
# Functions for setting model 0 MLE-based LRT adjustment parameters
param_func0 = lambda p0: np.concatenate([p0, [p0[1] * (1 - 1e-5), 0]])
step_func0 = lambda: params1[4:] * 0.01
bound_func0 = lambda p0: (np.array([0, 0]), np.array([p0[1], 1]))

param_func1 = lambda p0: np.concatenate([p0, [params1[4], 0]])
bound_func1 = lambda: (np.array([0, 0]), np.array([params0[1], 1]))

step_func2 = lambda: np.array([params1[5] * 0.01])
bounds2 = (np.array([0]), np.array([1]))


bound_func_model1 = lambda p1: (np.array([0, 0]), np.array([p1[1], 1]))


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


# Set up model arguments
_, __, model_args = dpluspy.uncerts.set_up_model_args(
    fit_file1, param_file1, bins=bins, pop_ids=pop_ids, u=u)


# Model 0 MLE
try:
    adj0_0 = dpluspy.uncerts.LRT_adjust(
        param_func0(params0), 
        model_args, 
        means, 
        varcovs, 
        bootreps, 
        nested_idx, 
        steps=step_func0(),
        bounds=bound_func0(params0),
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
    steps=step_func0(),
    bounds=bound_func1(),
    verbose=0,
)


adj0_2 = dpluspy.uncerts.LRT_adjust(
    param_func1(params0), 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx_2, 
    steps=step_func2(),
    bounds=bound_func1(),
    verbose=0,
)


# model 1 MLE
adj1 = dpluspy.uncerts.LRT_adjust(
    params1, 
    model_args, 
    means, 
    varcovs, 
    bootreps, 
    nested_idx, 
    bounds=bound_func_model1(params1),
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

