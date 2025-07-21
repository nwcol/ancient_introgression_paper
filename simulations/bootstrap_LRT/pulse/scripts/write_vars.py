
reps = list(range(100))
models = [0, 1]

with open("sim_vars.txt", "w") as fout:
    for model in models:
        for rep in reps:
            fout.write(f"{model},{rep}\n")