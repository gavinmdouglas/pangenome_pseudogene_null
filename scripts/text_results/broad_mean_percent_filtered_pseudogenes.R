rm(list = ls(all.names = TRUE))

# Mean and SD percent pseudogenes filtered out across all genomes.

filt_statsfiles <- list()

for (filename in list.files(path = '/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/pseudogenes/filt_intergenic_pseudo_stats',
                            pattern = '_stats.tsv',
                            full.names = TRUE)) {

  species <- gsub('_stats.tsv$', '', basename(filename))

  filt_statsfiles[[species]] <- read.table(filename, header = TRUE, sep = '\t', stringsAsFactors = FALSE)

}


for (filename in list.files(path = '/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/pseudogenes/filt_intergenic_pseudogene_stats',
                            pattern = '.txt',
                            full.names = TRUE)) {
  
  species <- gsub('.txt$', '', basename(filename))
  
  filt_statsfiles[[species]] <- read.table(filename, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
  
}

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

filt_statsfiles <- filt_statsfiles[names(filt_statsfiles)[which(names(filt_statsfiles) %in% rownames(pangenome))]]

filtered_stats <- do.call(rbind, filt_statsfiles)

filtered_stats$percent_failed <- ((filtered_stats$raw - filtered_stats$passed) / filtered_stats$raw) * 100

round(mean(filtered_stats$percent_failed), 2)
round(sd(filtered_stats$percent_failed), 2)
