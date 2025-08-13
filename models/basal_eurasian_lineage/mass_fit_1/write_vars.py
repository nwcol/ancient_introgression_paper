
nreps = 10
models = [
    "basal_eurasian_model_chag",
    "basal_eurasian_model_vin",
    "basal_eurasian_model_alt"
]
datasets = ["Bherer"]

with open("variables.txt", "w") as fout:
    for dataset in datasets:
        for model in models:
            for rep in range(nreps):
                fgraph = f"{model}.yaml"
                fparams = f"{model}_params.yaml"
                fdata = f"subset_{dataset}.pkl"
                fstem = f"{model}_{dataset}_{rep}"
                s = f"{fgraph},{fparams},{fdata},{fstem}\n"
                fout.write(s)
