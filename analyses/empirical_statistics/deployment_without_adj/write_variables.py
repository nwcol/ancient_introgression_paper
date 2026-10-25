
import pandas
import numpy as np

chroms = [str(x) for x in list(range(1, 23))]
df = pandas.read_csv(
    "../../../data/hg19.genome", names=["chrom", "L"], sep="\\s+")
Ldict = {c.replace("chr", ""): str(next(iter(df[df["chrom"] == c]["L"])))
         for c in df["chrom"]}

masks = [
    "accessible",
    "accessible_intergenic",
    "accessible_10kb_ccds_flank",
    "accessible_100kb_ccds_flank",
    "accessible_1e-4M_ccds_flank",
    "accessible_1e-3M_ccds_flank",
]

with open("variables.txt", "w") as fout:

    # Bherer maps
    for mask in masks:
        mask_path = "../../../data/bed_files/" + mask
        for chrom in chroms:
            args = [
                chrom,
                Ldict[chrom],
                mask_path,
                "../../../data/recombination_maps/Bherer",
                "sexavg_chr22.txt.gz",
                "pos",
                "cM",
                f"stats_{mask}_bherer_chr_{chrom}.pkl"
            ]
            fout.write(",".join(args) + "\n")

        for chrom in chroms:
            args = [
                chrom,
                Ldict[chrom],
                mask_path,
                "../../../data/recombination_maps/Zhou/JHS",
                "JHS.chr22.1000.plinkmap.txt.gz",
                "Position(bp)",
                "Map(cM)",
                f"stats_{mask}_zhou-JHS_chr_{chrom}.pkl"
            ]
            fout.write(",".join(args) + "\n")

    mask_subset = [masks[0], masks[1], masks[4]]
    for mask in mask_subset:
        mask_path = "../../../data/bed_files/" + mask
        for chrom in chroms:
            args = [
                chrom,
                Ldict[chrom],
                mask_path,
                "../../../data/recombination_maps/pyrho/YRI",
                f"YRI_recombination_map_hapmap_format_hg19_chr_22.txt.gz",
                "Position(bp)",
                "Map(cM)",
                f"stats_{mask}_pyrho-YRI_chr_{chrom}.pkl"
            ]
            fout.write(",".join(args) + "\n")

        for chrom in chroms:
            args = [
                chrom,
                Ldict[chrom],
                mask_path,
                "../../../data/recombination_maps/omni/YRI",
                f"YRI-{chrom}-final.txt.gz",
                "Position(bp)",
                "Map(cM)",
                f"stats_{mask}_omni-YRI_chr_{chrom}.pkl"
            ]
            fout.write(",".join(args) + "\n")



