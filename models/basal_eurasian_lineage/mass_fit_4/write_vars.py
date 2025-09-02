
nreps = 30
models = [
    "model_0_5",
    "model_0_6",
    "model_1_5",
    "model_1_6",
    "model_2_5",
    "model_2_6",
    "model_2_6_flipped",
    "model_3_5",
    "model_3_6",
]
datasets = ["Bherer", "ZhouJHS"]

with open("variables.txt", "w") as fout:
    for dataset in datasets:
        for model in models:
            for rep in range(nreps):
                args = [
                    f"{model}.yaml",
                    f"{model}_params.yaml",
                    f"subset_{dataset}.pkl",
                    f"{model}_{dataset}_fit_{rep}.yaml",
                ]
                varz = ",".join(args) + "\n"
                fout.write(varz)
