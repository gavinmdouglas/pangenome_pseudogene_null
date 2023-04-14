rm(list = ls(all.names = TRUE))

# Use  hard cut-off of <= 2 genomes for "ultra-rare" cloud with all other elements that fall in  <= 15% of genomes classified as "other cloud".
other.cloud_cutoff <- list()
soft.core_cutoff <- list()

pseudo_filt_genomes <- gsub(".txt$", "", list.files("/data1/gdouglas/projects/accessory_vs_pseudogene/pseudogenes/processed/pseudo_filtered_ids/"))

for (ACC_FILE in list.files(path = "/data1/gdouglas/projects/accessory_vs_pseudogene/genome_accessions/", full.names = FALSE)) {
  
  sp <- gsub("_accessions.txt", "", ACC_FILE)
  
  sp_genomes <- read.table(paste("/data1/gdouglas/projects/accessory_vs_pseudogene/genome_accessions/", ACC_FILE, sep = ""),
                           stringsAsFactors = FALSE)$V1
  
  num_filt_genomes <- length(which(sp_genomes %in% pseudo_filt_genomes))
  
  other.cloud_cutoff[[sp]] <- num_filt_genomes * 0.15
  soft.core_cutoff[[sp]] <- num_filt_genomes * 0.95
}

cluster_filt_lengths <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

# Restrict to clusters that remained after length filters.
cluster_breakdown <- cluster_breakdown[which(cluster_breakdown$cluster %in% rownames(cluster_filt_lengths)), ]

cluster_types <- readRDS(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_types.rds")

all_species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt",
                          header = FALSE, stringsAsFactors = FALSE)$V1


cluster_pangenome_categories <- list()

for (sp in all_species) {
  
  cluster_pangenome_categories[[sp]] <- list("pseudo" = list(),
                                             "intact" = list(),
                                             "both" = list())
  
  # Restrict to hits in specific species only and get table with just cluster and accession ids.
  sp_cluster_breakdown <- cluster_breakdown[which(cluster_breakdown$species == sp), ]
  sp_cluster_breakdown_simple <- sp_cluster_breakdown[, c("cluster", "accession")]
  
  # Remove lines representing multiple cluster instances in same genome.
  sp_cluster_breakdown_simple <- sp_cluster_breakdown_simple[-which(duplicated(sp_cluster_breakdown_simple)), ]
  
  cluster_frequency <- table(sp_cluster_breakdown_simple$cluster)
  
  cluster_frequency_intact <- cluster_frequency[which(names(cluster_frequency) %in% cluster_types$intact)]
  cluster_frequency_pseudo <- cluster_frequency[which(names(cluster_frequency) %in% cluster_types$intergenic.pseudogene)]
  cluster_frequency_both <- cluster_frequency[which(names(cluster_frequency) %in% cluster_types$both)]
  
  cluster_pangenome_categories[[sp]][["intact"]][["ultra.cloud"]] <- names(cluster_frequency_intact)[which(cluster_frequency_intact <= 2)]
  cluster_pangenome_categories[[sp]][["intact"]][["other.cloud"]] <- names(cluster_frequency_intact)[which(cluster_frequency_intact > 2 & cluster_frequency_intact <= other.cloud_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["intact"]][["soft.core"]] <- names(cluster_frequency_intact)[which(cluster_frequency_intact >= soft.core_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["intact"]][["shell"]] <- names(cluster_frequency_intact)[which(cluster_frequency_intact > other.cloud_cutoff[[sp]] & cluster_frequency_intact < soft.core_cutoff[[sp]])]
  
  cluster_pangenome_categories[[sp]][["pseudo"]][["ultra.cloud"]] <- names(cluster_frequency_pseudo)[which(cluster_frequency_pseudo <= 2)]
  cluster_pangenome_categories[[sp]][["pseudo"]][["other.cloud"]] <- names(cluster_frequency_pseudo)[which(cluster_frequency_pseudo > 2 & cluster_frequency_pseudo <= other.cloud_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["pseudo"]][["soft.core"]] <- names(cluster_frequency_pseudo)[which(cluster_frequency_pseudo >= soft.core_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["pseudo"]][["shell"]] <- names(cluster_frequency_pseudo)[which(cluster_frequency_pseudo > other.cloud_cutoff[[sp]] & cluster_frequency_pseudo < soft.core_cutoff[[sp]])]
  
  cluster_pangenome_categories[[sp]][["both"]][["ultra.cloud"]] <- names(cluster_frequency_both)[which(cluster_frequency_both <= 2)]
  cluster_pangenome_categories[[sp]][["both"]][["other.cloud"]] <- names(cluster_frequency_both)[which(cluster_frequency_both > 2 & cluster_frequency_both <= other.cloud_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["both"]][["soft.core"]] <- names(cluster_frequency_both)[which(cluster_frequency_both >= soft.core_cutoff[[sp]])]
  cluster_pangenome_categories[[sp]][["both"]][["shell"]] <- names(cluster_frequency_both)[which(cluster_frequency_both > other.cloud_cutoff[[sp]] & cluster_frequency_both < soft.core_cutoff[[sp]])]
  
}

# Save RDS containing breakdowns of which genes are in which categories.
saveRDS(object = cluster_pangenome_categories,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_pangenome_categories.rds")
