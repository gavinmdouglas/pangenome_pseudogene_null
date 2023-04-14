rm(list = ls(all.names = TRUE))

# Subset panaroo output tables to be the randomly selected 20 genomes for each species.

species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt",
                      stringsAsFactors = FALSE)$V1

for (sp in species) {

  subset_genomes_file <- outfile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/subsampled_genome_ids/',
                                          sp, '.txt', sep = '')

  subset_genomes <- read.table(file = subset_genomes_file, header = FALSE, stringsAsFactors = FALSE)$V1

  full_panaroo_file <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/focal_panaroo/',
                             sp, '/gene_presence_absence.csv.gz', sep = '')

  full_panaroo <- read.table(file = full_panaroo_file, header = TRUE, sep = ",", quote = "", stringsAsFactors = FALSE, check.names = FALSE)

  default_colnames <- colnames(full_panaroo)[1:3]
  
  subset_panaroo <- full_panaroo[, c(default_colnames, subset_genomes)]
  
  subset_panaroo <- subset_panaroo[which(rowSums(subset_panaroo[, subset_genomes] != '') > 0), ]
  
  subset_panaroo_outfile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/focal_panaroo/',
                                  sp, '/gene_presence_absence_ran.subset.csv', sep = '')
  write.table(x = subset_panaroo,
              file = subset_panaroo_outfile,
              quote = FALSE, sep = ',', row.names = FALSE, col.names = TRUE)
}
