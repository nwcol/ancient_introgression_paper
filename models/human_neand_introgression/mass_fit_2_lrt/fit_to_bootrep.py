
import sys
import dpluspy


graph_file =  sys.argv[1]
param_file = sys.argv[2]
data_file = sys.argv[3]
out_file = sys.argv[4]
rep_num = int(sys.argv[5])


# PARAMETERS
u = 1.3e-8

# Absolute difference in ll units required for cessation of model fitting
threshold = 1e-4
max_tries = 10
perturb = 0.0


pop_ids, bins, _, varcovs, bootreps = dpluspy.inference.load_bootstrap_reps(
    data_file, graph_file)
boot_means = bootreps[rep_num]


dpluspy.inference.optimize(
    graph_file,
    param_file,
    boot_means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    perturb=perturb,
    method="fmin",
    log=True,
    max_iter=2000,
    verbose=1,
    output=out_file,
    overwrite=True
)


# Fit the model to convergence using the Powell algorithm
converged = False
counter = 0
ll_log = []
while not converged:
    pnames, ps, ll = dpluspy.inference.optimize(
        out_file,
        param_file,
        boot_means,
        varcovs,
        pop_ids=pop_ids,
        bins=bins,
        u=u,
        method="powell",
        log=True,
        max_iter=100,
        verbose=500,
        output=out_file,
        overwrite=True
    )
    if counter > 0:
        last_ll = ll_log[-1]
        ll_diff = ll - last_ll
        if ll_diff < threshold:
            converged = True
        if counter > max_tries:
            converged = True
    ll_log.append(ll)
    counter += 1
