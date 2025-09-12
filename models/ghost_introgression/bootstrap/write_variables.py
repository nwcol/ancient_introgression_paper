
nreps = 100
datasets = ["Bherer", "ZhouJHS"]
models = [f"model_{i}" for i in range(4, 8)]

with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            for rep in range(nreps):
                vars = [
                    f"{model}_{dataset}_init.yaml",
                    f"{model}_params.yaml",
                    f"stats_{dataset}_weighted.pkl",
                    str(rep),
                    f"{model}_{dataset}_rep_{rep}.yaml"
                ]
                fout.write(",".join(vars) + "\n")