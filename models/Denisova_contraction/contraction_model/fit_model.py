
import dpluspy


model_file =  "model.yaml"
param_file = "params.yaml"
result0 = "fit_round0.yaml"
result1 = "fit_round1.yaml"


u = 1.3e-8
data_file = "../../../data/statistics/main/subset_stats.pkl"
pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph=model_file)

dpluspy.inference.optimize(
    model_file,
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    method="lbfgsb",
    log=True,
    max_iter=4000,
    verbose=10,
    output=result0,
)
dpluspy.inference.optimize(
    result0,
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
    output=result1,
)


model = dpluspy.inference.compute_bin_stats(result0, bins=bins, u=u, 
    sampled_demes=pop_ids)
dpluspy.plotting.plot_D_plus_curves(models=model, means=means, varcovs=varcovs,
    bins=bins, pop_ids=pop_ids, out="figure.pdf")