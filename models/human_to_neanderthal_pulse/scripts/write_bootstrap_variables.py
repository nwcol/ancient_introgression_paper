reps = list(range(100))
maps = ["Bherer", "ZhouFHS", "ZhouJHS"]
models = ["null_human_neanderthal", "human_neanderthal_pulse"]
with open("bootstrap_variables.txt", "w") as fout:
    for model in models:
        for mapp in maps:
            for rep in reps:
                s = (f"{model}_best_fit.yaml,{model}_params.yaml,"
                     f"subset_{mapp}.pkl,{rep},{model}_{mapp}_{rep}.csv\n")
                fout.write(s)