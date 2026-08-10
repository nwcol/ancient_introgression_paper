
nreps = 500
models = [
    "hn_model_1",
    "hn_model_5",
    "sd_model_1",
    "sd_model_5",
]


with open("sim_variables.txt", "w") as fout:
    for model in models:
        for rep in range(nreps):
            _args = [
                f"{model}.yaml", 
                f"{model}_rep_{rep}.pkl",
            ]
            args = ",".join(_args) + "\n"
            fout.write(args)