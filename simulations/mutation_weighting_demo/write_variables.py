with open("variables.txt", "w") as fout:
    for rep in range(5):
        fout.write(f"0,unif_rep_{rep}\n")
    for rep in range(5):
        fout.write(f"1,lmr_var_rep_{rep}\n")