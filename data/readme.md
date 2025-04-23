
The snakefile in this directory does not download all the needed input files itself, as some of them are very large. 
These, and their pathnames as expected by the snakefile, are listed here:

Simons genome project tarball
Four high-coverage archaic human genomes: "Vindija", "Chagyrskaya", "Altai" Neandertals, "Denisova" Denisovan
Three high-coverage ancient human genomes: "LBK", "Loschbour", "Ust'Ishim"


It builds the files required for empirical analysis using the D+ statistic, namely:

VCF files holding filtered/prepared genome sequences of seven ancient people and ---- modern African people
BED files corresponding to the above, designating callable areas of the genome and used to compute the denominators of statistics
Recombination map files
Site-resolution mutation maps (Roulette)
Files representing genomic windows for the block-bootstrap


