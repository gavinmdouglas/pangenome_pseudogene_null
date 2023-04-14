rm(list = ls(all.names = TRUE))

nonfocal <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/metadata_highqual_RefSeq_non.focal_max20.rds")

intersect(grep("Bacteria", nonfocal$gtdb_taxonomy), grep("Archaea", nonfocal$gtdb_taxonomy))

bacteria_nonfocal <- nonfocal[grep("Bacteria", nonfocal$gtdb_taxonomy), ]
archaea_nonfocal <- nonfocal[grep("Archaea", nonfocal$gtdb_taxonomy), ]

bacteria_nonfocal_nonfragmented_ids <- rownames(bacteria_nonfocal[which(bacteria_nonfocal$contig_count < 50), ])
bacteria_nonfocal_fragmented_ids <- rownames(bacteria_nonfocal[which(bacteria_nonfocal$contig_count >= 50), ])

write.table(x = bacteria_nonfocal_nonfragmented_ids,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/bacteria_nonfragmented_ids.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

write.table(x = bacteria_nonfocal_fragmented_ids,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/bacteria_genomes_fragmented_ids.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")


archaea_nonfocal_nonfragmented_ids <- rownames(archaea_nonfocal[which(archaea_nonfocal$contig_count < 50), ])
archaea_nonfocal_fragmented_ids <- rownames(archaea_nonfocal[which(archaea_nonfocal$contig_count >= 50), ])

write.table(x = archaea_nonfocal_nonfragmented_ids,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/archaea_nonfragmented_ids.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")

write.table(x = archaea_nonfocal_fragmented_ids,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/archaea_genomes_fragmented_ids.txt",
            quote = FALSE, col.names = FALSE, row.names = FALSE, sep = "\t")
