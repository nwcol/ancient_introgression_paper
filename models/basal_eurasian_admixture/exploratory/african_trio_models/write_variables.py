datasets = ["Bherer"]
models = ["model_0", "model_1", "model_2"]
with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            for rep in range(5):
                fgraph = f"{model}.yaml"
                fparams = f"{model}_params.yaml"
                fdata = f"subset_{dataset}.pkl"
                foutstem = f"{model}_rep{rep}"
                s = f"{fgraph},{fparams},{fdata},{foutstem}\n"
                fout.write(s)