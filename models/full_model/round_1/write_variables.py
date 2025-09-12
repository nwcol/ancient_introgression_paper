
nreps = 50
datasets = ["Bherer", "ZhouJHS"]

with open("variables.txt", "w") as fout:
    for dataset in datasets:
        for rep in range(nreps):
            vars = [
                f"full_model_round_0_init.yaml",
                f"full_model_round_0_params.yaml",
                f"stats_{dataset}_weighted.pkl",
                f"full_model_round_0_{dataset}_fit_{rep}.yaml"
            ]
            fout.write(",".join(vars) + "\n")