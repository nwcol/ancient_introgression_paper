
import pandas

df = pandas.read_csv("../../../data/hg19.genome", names=["chrom", "L"], sep="\\s+")
Ldict = {int(c.replace("chr", "")): next(iter(df[df["chrom"] == c]["L"]))
         for c in df["chrom"]}

reps = 10
with open("variables.txt", "w") as fout:
    for i in range(reps):
        for c in range(1,23):
            line = ",".join(
                [str(c),
                str(Ldict[c]),
                "graph.yaml",
                "pops.txt",
                f"rep_{i}_chr_{c}"]) + "\n"
            fout.write(line)

