
nreps = 30
datasets = ["Bherer"]
models = [f"model_{i}" for i in range(8)]

with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            for rep in range(nreps):
                vars = [
                    f"{model}_init.yaml",
                    f"{model}_params.yaml",
                    f"subset_{dataset}.pkl",
                    f"{model}_fit_{rep}.yaml"
                ]
                fout.write(",".join(vars) + "\n")