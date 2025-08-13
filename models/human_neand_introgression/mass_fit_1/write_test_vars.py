
nreps = 1
datasets = ["Bherer"]
models = [
    "human_neand_model_0",
    "human_neand_model_1",
    "human_neand_model_2",
    "human_neand_model_3",
    "human_neand_model_4"
]
samples = ["yor1", "mbu1"]

with open("test_variables.txt", "w") as fout:
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