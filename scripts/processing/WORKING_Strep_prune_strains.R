rm(list = ls(all.names = TRUE))

# Prune Strep pneu strains that are most similar to get a set of 1000 genomes.
# Needed so that panaroo is run on a more feasible dataset.

Strep_ani <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/dashing/Streptococcus_pneumoniae.output.tsv.gz",
                        header = TRUE, row.names = 1, sep = "\t", comment.char = "", na.strings = "-")

rownames(Strep_ani) <- gsub(".fna$", "", rownames(Strep_ani))
colnames(Strep_ani) <- gsub(".fna$", "", colnames(Strep_ani))

rownames(Strep_ani) <- gsub("^/data1/gdouglas/projects/accessory_vs_pseudogene/genomes/Streptococcus_pneumoniae/", "", rownames(Strep_ani))
colnames(Strep_ani) <- gsub("^X.data1.gdouglas.projects.accessory_vs_pseudogene.genomes.Streptococcus_pneumoniae.", "", colnames(Strep_ani))