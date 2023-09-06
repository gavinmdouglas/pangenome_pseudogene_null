rm(list = ls(all.names = TRUE))

# Get rsync download commands for species outside of the 10 focal ones.

# Restrict to species with at least 10 genomes, and take a maximum of 20.

bacteria_metadata <- read.table("/data1/gdouglas/db/gtdbtk_data_r202/bac120_metadata_r202.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

archaea_metadata <- read.table("/data1/gdouglas/db/gtdbtk_data_r202/ar122_metadata_r202.tsv.gz",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

combined_metadata <- rbind(archaea_metadata, bacteria_metadata)

rownames(combined_metadata) <- gsub("^GB_", "", rownames(combined_metadata))
rownames(combined_metadata) <- gsub("^RS_", "", rownames(combined_metadata))

high.qual_ids <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/highqual_gtdb_genomes.txt",
                            header = FALSE, stringsAsFactors = FALSE)$V1

combined_metadata_highqual <- combined_metadata[high.qual_ids, ]

# Subset to genomes in RefSeq only.
combined_metadata_highqual_RefSeq <- combined_metadata_highqual[grep("^GCF_", rownames(combined_metadata_highqual)), ]



# Make sure that the 10 focal species are ignored.
focal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt",
                            stringsAsFactors = FALSE)$V1

focal_species_accessions <- read.table(file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/genome_info/accessions.tsv.gz',
                                       stringsAsFactors = FALSE, header = TRUE, sep = '\t')$V2


combined_metadata_highqual_RefSeq_non.focal <- combined_metadata_highqual_RefSeq[which(! rownames(combined_metadata_highqual_RefSeq) %in% focal_species_accessions), ]

combined_metadata_highqual_RefSeq_non.focal

species_tallies <- table(combined_metadata_highqual_RefSeq_non.focal$gtdb_taxonomy)
length(species_tallies)

species_tallies_min10 <- species_tallies[which(species_tallies >= 10)]
length(species_tallies_min10)

species_tallies_min10_max20 <- species_tallies_min10
species_tallies_min10_max20[which(species_tallies_min10_max20 > 20)] <- 20
length(species_tallies_min10_max20)

# So there should be 660 non-focal additional species to download (after randomly subsetting them all to 20 genomes if needed).

# Subsample table.
species_20_or_fewer_genomes <- names(species_tallies_min10)[which(species_tallies_min10 <= 20)]
species_greater_than_20_genomes <- names(species_tallies_min10)[which(species_tallies_min10 > 20)]

combined_metadata_highqual_RefSeq_non.focal_below20 <- combined_metadata_highqual_RefSeq_non.focal[which(combined_metadata_highqual_RefSeq_non.focal$gtdb_taxonomy %in% species_20_or_fewer_genomes), ]

combined_metadata_highqual_RefSeq_non.focal_above20_RAW <- list()

for (tmp_species in species_greater_than_20_genomes) {
  tmp_species_subset <- combined_metadata_highqual_RefSeq_non.focal[which(combined_metadata_highqual_RefSeq_non.focal$gtdb_taxonomy == tmp_species), ]
  combined_metadata_highqual_RefSeq_non.focal_above20_RAW[[tmp_species]] <- tmp_species_subset[sample(x = 1:nrow(tmp_species_subset), size = 20), ]
}

combined_metadata_highqual_RefSeq_non.focal_above20 <- do.call(rbind, combined_metadata_highqual_RefSeq_non.focal_above20_RAW)
rownames(combined_metadata_highqual_RefSeq_non.focal_above20) <- gsub("^.*\\.GCF_", "GCF_", rownames(combined_metadata_highqual_RefSeq_non.focal_above20))


combined_metadata_highqual_RefSeq_non.focal_max20 <- rbind(combined_metadata_highqual_RefSeq_non.focal_above20, combined_metadata_highqual_RefSeq_non.focal_below20)

# Get mapping of ids to use.
combined_metadata_highqual_RefSeq_non.focal_max20$file_prefix <- paste(rownames(combined_metadata_highqual_RefSeq_non.focal_max20), combined_metadata_highqual_RefSeq_non.focal_max20$ncbi_assembly_name, sep = "_")
combined_metadata_highqual_RefSeq_non.focal_max20$file_prefix <- gsub(" ", "_", combined_metadata_highqual_RefSeq_non.focal_max20$file_prefix)

# Tip ids input to tanger cant have underscores or periods, so make a mapping for these tips too while I'm at it.
combined_metadata_highqual_RefSeq_non.focal_max20$ranger_tip_id <- gsub("_", "", rownames(combined_metadata_highqual_RefSeq_non.focal_max20))
combined_metadata_highqual_RefSeq_non.focal_max20$ranger_tip_id <- gsub("\\.", "", combined_metadata_highqual_RefSeq_non.focal_max20$ranger_tip_id)

combined_metadata_highqual_RefSeq_non.focal_max20$accession <- rownames(combined_metadata_highqual_RefSeq_non.focal_max20)

write.table(x = combined_metadata_highqual_RefSeq_non.focal_max20[, c("accession", "file_prefix", "ranger_tip_id")],
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/non.focal.subsampled_accession_id_map.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")


# Also generate all the rsync download commands that are needed.

# E.g.
# rsync --copy-links --times  rsync://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/902/209/435/GCF_902209435.1_Katie-F-Assembly-1/GCF_902209435.1_Katie-F-Assembly-1_genomic.fna.gz ./


command_start = "rsync --copy-links --times rsync://ftp.ncbi.nlm.nih.gov/genomes/all"

rsync_cmds <- c()

for (accession in rownames(combined_metadata_highqual_RefSeq_non.focal_max20)) {
  accession_split <- strsplit(accession, "_")[[1]]
  
  accession_id_split <- strsplit(accession_split[2], "")[[1]]
  
  id_part1 <- paste0(accession_id_split[1:3], collapse = "")
  id_part2 <- paste0(accession_id_split[4:6], collapse = "")
  id_part3 <- paste0(accession_id_split[7:9], collapse = "")
  
  file_prefix <- combined_metadata_highqual_RefSeq_non.focal_max20[accession, "file_prefix"]
  
  fna_cmd_ending <- paste0(file_prefix, "_genomic.fna.gz ./")
  md5_cmd_ending <- paste0("md5checksums.txt ", "./", file_prefix, "_genomic.fna.gz.exp.md5")
  
  fna_cmd <- paste(command_start, accession_split[1], id_part1, id_part2, id_part3, file_prefix, fna_cmd_ending, sep = "/")
  md5_cmd <- paste(command_start, accession_split[1], id_part1, id_part2, id_part3, file_prefix, md5_cmd_ending, sep = "/")
  
  rsync_cmds <- c(rsync_cmds, fna_cmd, md5_cmd)
}

write.table(x = rsync_cmds,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/genome_rsync_download_cmds.sh",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

saveRDS(object = combined_metadata_highqual_RefSeq_non.focal_max20, file = "/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/metadata_highqual_RefSeq_non.focal_max20.rds")
