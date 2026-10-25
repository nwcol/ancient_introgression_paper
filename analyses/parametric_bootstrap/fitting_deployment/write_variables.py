
nreps = 500
ntrials = 5
with open("variables.txt", "w") as fout:
    for i in range(nreps):
        for j in range(ntrials):
            fout.write(f"{i},{j}\n")
