
nreps = 100
chroms = list(range(1, 23))

with open("sim_variables", "w") as fout:
    for rep in range(nreps):
        for chrom in chroms:
            _args = [
                str(chrom),
                f"scratch_chr_{chrom}.vcf",
                f"control_rep_{rep}_chr_{chrom}_stats.pkl",
            ]
            args = ",".join(_args) + "\n"
            fout.write(args)
