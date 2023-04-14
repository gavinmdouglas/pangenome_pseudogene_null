rm(list = ls(all.names = TRUE))

nonfocal <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/metadata_highqual_RefSeq_non.focal_max20.rds")

nonfocal$species <- gsub("^.*;s__", "", nonfocal$gtdb_taxonomy)

for (sp in unique(nonfocal$species)) {
  sp_ids <- rownames(nonfocal)[which(nonfocal$species == sp)]
  sp_nospace <- gsub(" ", "_", sp)
  write.table(x = sp_ids,
              file = paste("/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/ids_by_species/", sp_nospace, ".txt", sep = ""),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}



