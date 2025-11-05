
chroms = list(range(1, 23))
masks = [
    "intergenic",
    "10kb_exon_buffer",
    "20kb_exon_buffer",
    "50kb_exon_buffer",
    "1e-4M_exon_buffer"
]
with open("buffer_variables.txt", "w") as fout:
    for mask in masks:
        for chrom in chroms:
            fout.write(f"{mask},{chrom}\n")
