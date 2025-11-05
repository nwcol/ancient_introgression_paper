
with open("variables.txt", "w") as fout:
    for i in range(4):
        for j in range(1, 23):
            fout.write(f"{i},{j}\n")