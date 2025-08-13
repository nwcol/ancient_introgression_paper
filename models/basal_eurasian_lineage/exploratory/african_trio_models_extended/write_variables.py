datasets = ["Bherer"]
models = ["model_ust0", "model_los0", "model_stu0",
          "model_ust1", "model_los1", "model_stu1"]
with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            fgraph = f"{model}.yaml"
            fparams = f"{model}_params.yaml"
            fdata = f"subset_{dataset}.pkl"
            foutstem = f"{model}_fit"
            s = f"{fgraph},{fparams},{fdata},{foutstem}\n"
            fout.write(s)