
models = [
    ("equilibrium.yaml", "pops1.txt"),
    ("growth_2.yaml", "pops1.txt"),
    ("growth_50.yaml", "pops1.txt"),
    ("shrink_0.2.yaml", "pops1.txt"),
    ("shrink_5.yaml", "pops1.txt"),
    ("isolation_0.yaml", "pops2.txt"),
    ("isolation_1.yaml", "pops2.txt"),
    ("isolation_4.yaml", "pops2.txt"),
    ("admixture_2_0.yaml", "pops2.txt"),
    ("admixture_2_0.05.yaml", "pops2.txt"),
    ("admixture_4_0.yaml", "pops2.txt"),
    ("admixture_4_0.05.yaml", "pops2.txt"),
]

with open("variables.txt", "w") as fout:
    for model, pops in models:
        for rep in range(200):
            outname = model.replace(".yaml", "") + f"_rep_{rep}.pkl"
            vars = [
                model,
                pops,
                outname
            ]
            line = ",".join(vars) + "\n"
            fout.write(line)

