
nreps = 30
datasets = ["Bherer", "ZhouJHS"]
models = [f"model_{i}" for i in range(3)]

with open("variables.txt", "w") as fout:
    for model in models:
        for dataset in datasets:
            for rep in range(nreps):
                vars = [
                    f"{model}_init.yaml",
                    f"{model}_params.yaml",
                    f"stats_{dataset}_weighted.pkl",
                    f"{model}_{dataset}_fit_{rep}.yaml"
                ]
                fout.write(",".join(vars) + "\n")