reps = list(range(5))
maps = ["Bherer", "ZhouFHS", "ZhouJHS"]
models = ["null_human_neanderthal", "human_neanderthal_pulse"]
with open("variables.txt", "w") as fout:
    for model in models:
        for mapp in maps:
            for rep in reps:
                s = (f"{model}.yaml,{model}_params.yaml,subset_{mapp}.pkl,"
                    f"{model}_{mapp}_{rep}\n")
                fout.write(s)