
reps = list(range(100))
data_sets = ["Bherer", "ZhouJHS"]
models = ["model_0", "model_1", "model_2"]
with open("lrt_variables.txt", "w") as fout:
    for model in models:
        for data_set in data_sets:
            for rep in reps:
                graph_file = f"{model}_init.yaml"
                param_file = f"{model}_params.yaml"
                data_file = f"subset_{data_set}.pkl"
                out_file = f"{model}_{data_set}_rep_{rep}.yaml"
                s = ",".join([graph_file, param_file, data_file, out_file, str(rep)]) + "\n"
                fout.write(s)