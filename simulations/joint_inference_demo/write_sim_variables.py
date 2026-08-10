with open("sim_variables.txt", "w") as fout:
    for rep in range(100, 500):
        fout.write(f"rep_{rep}.pkl\n")
