
nreps = 500

with open("fit_variables.txt", "w") as fout:
    for model in ["model_1", "model_5"]:
        for rep in range(nreps):
            vars = [
                "hn_null_model.yaml",
                "hn_null_model_params.yaml",
                "hn_alt_model.yaml",
                "hn_alt_model_params.yaml",
                f"hn_{model}_rep_{rep}.pkl",
                f"hn_{model}_rep_{rep}_fit"
            ]
            fout.write(",".join(vars) + "\n")

    for model in ["model_1", "model_5"]:
        for rep in range(nreps):
            vars = [
                "sd_null_model.yaml",
                "sd_null_model_params.yaml",
                "sd_alt_model.yaml",
                "sd_alt_model_params.yaml",
                f"sd_{model}_rep_{rep}.pkl",
                f"sd_{model}_rep_{rep}_fit"
            ]
            fout.write(",".join(vars) + "\n")