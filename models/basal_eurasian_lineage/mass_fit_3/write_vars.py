
nreps = 30
models = [
    "null_model",
    "admix_model",
    "deep_struc_model",
    "shallow_struc_model"
]
datasets = ["Bherer"]

with open("variables.txt", "w") as fout:
    for dataset in datasets:
        for model in models:
            for rep in range(nreps):
                args = [
                    f"{model}.yaml",
                    f"{model}_params.yaml",
                    f"subset_{dataset}.pkl",
                    f"{model}_fit_{rep}.yaml",
                ]
                varz = ",".join(args) + "\n"
                fout.write(varz)
