
nreps = 10
models = [
    "null_model_chag",
    "null_model_vin",
    "admix_model_chag",
    "admix_model_vin",
    "marg_model_chag",
    "marg_model_vin"
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
