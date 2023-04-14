rm(list = ls(all.names = TRUE))

bacteria_metadata <- read.table("/data1/gdouglas/db/gtdbtk_data_r202/bac120_metadata_r202.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

archaea_metadata <- read.table("/data1/gdouglas/db/gtdbtk_data_r202/ar122_metadata_r202.tsv.gz",
                               header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

combined_metadata <- rbind(archaea_metadata, bacteria_metadata)

rownames(combined_metadata) <- gsub("^GB_", "", rownames(combined_metadata))
rownames(combined_metadata) <- gsub("^RS_", "", rownames(combined_metadata))

combined_metadata_highmimag <- combined_metadata[which(combined_metadata$mimag_high_quality == "t"), ]
combined_metadata_extra_highqual <- combined_metadata_highmimag[which(combined_metadata_highmimag$checkm_completeness > 98 & combined_metadata_highmimag$checkm_contamination < 1), ]
combined_metadata_extra_highqual <- combined_metadata_extra_highqual[which(combined_metadata_extra_highqual$contig_count < 1000), ]
combined_metadata_extra_highqual <- combined_metadata_extra_highqual[which(combined_metadata_extra_highqual$n50_contigs > 5000), ]
combined_metadata_extra_highqual <- combined_metadata_extra_highqual[which(combined_metadata_extra_highqual$ambiguous_bases < 100000), ]

# Needed to also remove GCF_000702925.2 as this was suppressed in RefSeq/GenBank
combined_metadata_extra_highqual <- combined_metadata_extra_highqual[-which(rownames(combined_metadata_extra_highqual) == "GCF_000702925.2"), ]

combined_metadata_extra_highqual$sp <- gsub('^.*;s__', '', combined_metadata_extra_highqual$gtdb_taxonomy)
combined_metadata_extra_highqual$sp <- gsub(' ', '_', combined_metadata_extra_highqual$sp)


# Subset to table to species used:
focal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/focal_species.txt",
                            stringsAsFactors = FALSE, header = FALSE)$V1
nonfocal_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/non.focal_species.txt",
                               stringsAsFactors = FALSE, header = FALSE)$V1
focal_and_nonfocal_species <- sort(c(focal_species, nonfocal_species))

taxa_working <- combined_metadata_extra_highqual[which(combined_metadata_extra_highqual$sp %in% focal_and_nonfocal_species), c("sp", "gtdb_taxonomy")]
taxa_working <- taxa_working[-which(duplicated(taxa_working)), ]
dim(taxa_working)

taxa_working$phylum <- gsub(';c__.*$', '', gsub('^.*;p__', 'p__', taxa_working$gtdb_taxonomy))
taxa_working$class <- gsub(';o__.*$', '', gsub('^.*;c__', 'c__', taxa_working$gtdb_taxonomy))
taxa_working$order <- gsub(';f__.*$', '', gsub('^.*;o__', 'o__', taxa_working$gtdb_taxonomy))
taxa_working$family <- gsub(';g__.*$', '', gsub('^.*;f__', 'f__', taxa_working$gtdb_taxonomy))
taxa_working$genus <- gsub(';s__.*$', '', gsub('^.*;g__', 'g__', taxa_working$gtdb_taxonomy))

rownames(taxa_working) <- NULL

write.table(x = taxa_working,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/focal_and_non.focal_taxonomy.tsv",
            quote = FALSE, col.names = TRUE, row.names = FALSE, sep = '\t')

