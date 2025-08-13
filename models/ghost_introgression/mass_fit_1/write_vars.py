
nreps = 30
datasets = ["Bherer", "ZhouJHS"]
models = [f"ghost_model_{i}" for i in range(16)]
samples = ["yor1"]

with open("variables.txt", "w") as fout:
    for model in models:
        for sample in samples:
            for dataset in datasets:
                for rep in range(nreps):
                    vars = [
                        f"{model}_{sample}.yaml",
                        f"{model}_{sample}_params.yaml",
                        f"subset_{dataset}.pkl",
                        f"{model}_{dataset}_{sample}_rep_{rep}.yaml"
                    ]
                    fout.write(",".join(vars) + "\n")