# Write ghost-Denisova introgression grid models.

import numpy as np
import demes 
import moments


param_file = "model_builders/ghost_params.yaml"
options = moments.Demes.Inference._get_params_dict(param_file)

model1_file = "model_builders/model1.yaml"
model2_file = "model_builders/model2.yaml"


T0s = [9e5, 1e6, 1.1e6, 1.3e6, 1.5e6, 1.7e6, 1.9e6, 2.2e6, 2.5e6]
T1s = [1e5, 1.5e5, 2e5, 2.5e5, 3e5, 4e5, 5e5]
ps = [0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10]


def write_model(params, in_file, out_file):
    builder = moments.Demes.Inference._get_demes_dict(in_file)
    builder = moments.Demes.Inference._update_builder(builder, options, params)
    graph = demes.Graph.fromdict(builder)
    demes.dump(graph, out_file)


counter = 1 
for model in (model1_file, model2_file):
    for T0 in T0s:
        for T1 in T1s:
            for pp in ps:
                params = np.array([T0, T1, pp])
                out_file = f"models/grid_model_{counter}.yaml"
                write_model(params, model, out_file)
                counter += 1



