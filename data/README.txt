The dataset used in this project can be replicated by running Snakefiles in
the subdirectories of this directory. Some very large files are not downloaded
by the Snakefiles- you must download these manually (e.g. using wget), place 
them where you wish, and edit config.yaml so that it contains the correct file
paths. These large files are:

- SGDP tarball
- contents of the following directories from the Max Planck ftp:
    http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Altai/
    http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Denisova/
    http://ftp.eva.mpg.de/neandertal/Vindija/VCF/LBK/
    http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Loschbour/
    http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Ust_Ishim/
    http://ftp.eva.mpg.de/neandertal/Vindija/VCF/Vindija33.19/
    http://ftp.eva.mpg.de/neandertal/Chagyrskaya/VCF/
- the Roulette recombination map in VCF form, build GRCh37:
    http://genetics.bwh.harvard.edu/downloads/Vova/Roulette/hg19/

With these files secured, you will want to run these Snakefiles sequentially:
- mutation_maps/Snakefile 
- bed_files/Snakefile
- intervals/Snakefile [TODO write this one!]
These can be executed in any order, before or after the above:
- recombination_maps/Snakefile 
- vcf_files/Snakefile 
With all that done, run the Snakefile in statistics/ to compute the D+ statistic
from this data under the configurations used in this work.

In places I have been a little sloppy with directory and file names. Note that 
throughout I use several synonyms for GRCh37!


