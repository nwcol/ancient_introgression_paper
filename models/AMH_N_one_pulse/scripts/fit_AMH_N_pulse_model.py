
import sys
import dpluspy


trial_id = sys.argv[1]

graph_file = "pulse.model.yaml"
param_file = "pulse.params.yaml"
output_0 = f"pulse.fitted.{trial_id}_rep0.yaml"
output_1 = f"pulse.fitted.{trial_id}_rep1.yaml"
u = 1.3e-8
data_file = "subset_stats.pkl"


pop_ids, bins, means, varcovs = dpluspy.inference.load_stats(
    data_file, graph_file)

# Round 0 of inference, using 5000x fmin and large perturbation
dpluspy.inference.optimize(
    graph_file,
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    perturb=0.5,
    method="fmin",
    max_iter=5000,
    verbose=50,
    output=output_0
)
# Round 2 using 100x powell
dpluspy.inference.optimize(
    output_0,
    param_file,
    means,
    varcovs,
    pop_ids=pop_ids,
    bins=bins,
    u=u,
    method="powell",
    max_iter=100,
    verbose=20,
    output=output_1
)