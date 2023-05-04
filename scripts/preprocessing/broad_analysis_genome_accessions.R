rm(list = ls(all.names = TRUE))

# Create file of all genome accessions used for the broad pangenome analysis.
# Not all of these genomes were able to be downloaded, so include a column
# indicating whether it was downloaded or not.

nonfocal <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/metadata_highqual_RefSeq_non.focal_max20.rds")

accessions <- data.frame(species = gsub(' ', '_', gsub("^.*;s__", "", nonfocal$gtdb_taxonomy)),
                         accession = rownames(nonfocal))

nonfocal_annotated_accessions <- gsub('.gff', '', list.files("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/non.focal_genomes/non.focal_gffs/"))

accessions$could_download <- FALSE
accessions$could_download[which(accessions$accession %in% nonfocal_annotated_accessions)] <- TRUE

focal_accessions_raw <- list()
for (focal_accession_file in list.files("/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/subsampled_genome_ids/", full.names = TRUE)) {
  species <- gsub('.txt$', '', basename(focal_accession_file))
  species_accessions <- read.table(focal_accession_file, stringsAsFactors = FALSE)$V1
  focal_accessions_raw[[species]] <- data.frame(species = species,
                                                accession = species_accessions,
                                                could_download = TRUE)
}

focal_accessions_df <- do.call(rbind, focal_accessions_raw)
rownames(focal_accessions_df) <- NULL

accessions <- rbind(accessions, focal_accessions_df)

write.table(x = accessions,
            file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/accessions.tsv",
            row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
