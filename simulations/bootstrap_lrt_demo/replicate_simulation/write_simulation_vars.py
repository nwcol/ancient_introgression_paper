
nreps = 10
nsamps = 100
models = [
    "hn_null_model",
    "hn_model_1",
    "hn_model_5",
    "sd_null_model",
    "sd_model_1",
    "sd_model_5",
]


with open("sim_vars.txt", "w") as fout:
    for model in models:
        for rep in range(nreps):
            for sample in range(nsamps):
                _args = [
                    f"{model}_rep_{rep}.yaml", 
                    f"{model}_rep_{rep}_samp_{sample}.pkl",
                ]
                args = ",".join(_args) + "\n"
                fout.write(args)