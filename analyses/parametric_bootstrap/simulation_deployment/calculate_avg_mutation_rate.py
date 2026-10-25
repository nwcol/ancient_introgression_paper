
import numpy as np
import h2py


input_maps = [f"data/map_chr{c}.tsv.gz" for c in range(1, 23)]
input_masks = [f"data/mask_chr{c}.bed.gz" for c in range(1, 23)]

n_sites = 0
tot_u = 0.0

for map_file, mask_file in zip(input_maps, input_masks):
    mask = h2py.GeneticMask.from_bed_file(mask_file)
    positions = mask.to_positions()
    umap = h2py.parsing._get_mut_rates(map_file, positions)
    tot_u += np.sum(umap)
    n_sites += len(positions)
    print(f"{map_file}, {mask_file}: {n_sites} sites")

print(f"avg mutation rate: {tot_u / n_sites:.6}")
