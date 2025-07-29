datasets = ["Bherer"]
models = ["model_ust", "model_los", "model_stu"]
with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            fgraph = f"{model}.yaml"
            fparams = f"{model}_params.yaml"
            fdata = f"subset_{dataset}.pkl"
            foutstem = f"{model}_fit"
            s = f"{fgraph},{fparams},{fdata},{foutstem}\n"
            fout.write(s)