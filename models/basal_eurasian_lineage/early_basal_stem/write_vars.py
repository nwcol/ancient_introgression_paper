
nreps = 50
models = [
    "basal_stem_2k",
    "basal_stem_5k",
    "basal_stem_10k",
    "basal_stem_S",
    "basal_stem_MH"
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
