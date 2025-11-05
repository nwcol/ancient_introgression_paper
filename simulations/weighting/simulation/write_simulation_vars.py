
nreps = 100


with open("test_vars.txt", "w") as fout:
    for heterog_u in [0, 1]:
        for rep in range(1):
            _args = [
                "full_model.yaml",
                f"heterog_u_{heterog_u}_rep_{rep}",
                str(heterog_u),
            ]
            args = ",".join(_args) + "\n"
            fout.write(args)


with open("sim_vars.txt", "w") as fout:
    for heterog_u in [0, 1]:
        for rep in range(nreps):
            _args = [
                "full_model.yaml",
                f"heterog_u_{heterog_u}_rep_{rep}",
                str(heterog_u),
            ]
            args = ",".join(_args) + "\n"
            fout.write(args)