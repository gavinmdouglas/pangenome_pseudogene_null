rm(list = ls(all.names = TRUE))

# Identify sets of high quality genomes to download.

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

write.table(x = rownames(combined_metadata_extra_highqual),
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/highqual_gtdb_genomes.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE)
