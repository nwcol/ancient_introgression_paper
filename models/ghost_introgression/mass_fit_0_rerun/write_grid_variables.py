

with open("grid_variables.txt", "w") as fout:
    for ii in range(693):
        model = f"ghost_den_start_{ii}.yaml"
        output = f"ghost_den_fit_{ii}"
        s = f"{model},{output}\n"
        fout.write(s)