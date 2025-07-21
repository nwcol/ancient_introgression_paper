

# These are the filename pattern, pos. col name, map col name
map_types = {
    "Bherer": ("sexavg_chr{chrom}.txt.gz", "pos", "cM"),
    "hapmap": ("genetic_map_GRCh37_chr{chrom}.txt.gz", "Position(bp)", "Map(cM)"),
    "Hinch": ("maps_chr.{chrom}.txt.gz", "Physical_Pos", "AA_Map"),
    "omniFIN": ("FIN-{chrom}-final.txt.gz", "Position(bp)", "Map(cM)"),
    "omniLWK": ("LWK-{chrom}-final.txt.gz", "Position(bp)", "Map(cM)"),
    "omniYRI": ("YRI-{chrom}-final.txt.gz", "Position(bp)", "Map(cM)"),
    "pyrhoFIN": ("FIN_recombination_map_hapmap_format_hg19_chr_{chrom}.txt.gz", "Position(bp)", "Map(cM)"),
    "pyrhoLWK": ("LWK_recombination_map_hapmap_format_hg19_chr_{chrom}.txt.gz", "Position(bp)", "Map(cM)"),
    "pyrhoYRI": ("YRI_recombination_map_hapmap_format_hg19_chr_{chrom}.txt.gz", "Position(bp)", "Map(cM)"),
    "ZhouFHS": ("FHS.chr{chrom}.1000.plinkmap.txt.gz", "Position(bp)", "Map(cM)"),
    "ZhouJHS": ("JHS.chr{chrom}.1000.plinkmap.txt.gz", "Position(bp)", "Map(cM)") 
}


with open("variables", "w") as fout:
    for name in map_types:
        pattern, pos_col, map_col = map_types[name]
        for chrom in range(1, 23):
            filename = pattern.format(chrom=chrom)
            line = ",".join([str(chrom), name, filename, pos_col, map_col])+"\n"
            fout.write(line)

