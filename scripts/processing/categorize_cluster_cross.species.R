rm(list = ls(all.names = TRUE))

# Get breakdown of what clusters are shared across different species (and what segment of each species' pangenome it is found in).
# May not actually include this in the manuscript, but good to have in case.

cluster_filt_lengths <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

# Restrict to clusters that remained after length filters.
cluster_breakdown <- cluster_breakdown[which(cluster_breakdown$cluster %in% rownames(cluster_filt_lengths)), ]

cluster_pangenome_categories <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_pangenome_categories.rds")


# Identify cross-species clusters and summarize what species they are in and what category subtype they are (e.g., "both_cloud").
cluster_breakdown_by_species <- cluster_breakdown[, c("species", "cluster")]

cluster_breakdown_single_by_species <- cluster_breakdown[-which(duplicated(cluster_breakdown_by_species)), ]

cluster_breakdown_single_by_species_tally <- table(cluster_breakdown_single_by_species$cluster)

cross_species_clusters <- names(cluster_breakdown_single_by_species_tally)[which(cluster_breakdown_single_by_species_tally > 1)]

cluster_breakdown_single_by_species_cross.species <- cluster_breakdown_single_by_species[which(cluster_breakdown_single_by_species$cluster %in% cross_species_clusters), ]


cross_species_cluster_summary <- list()

for (cross_species_cluster in cross_species_clusters) {
  
  cross_species_cluster_subset <- cluster_breakdown_single_by_species_cross.species[which(cluster_breakdown_single_by_species_cross.species$cluster == cross_species_cluster), ]
  
  category_vec <- as.character()
  
  for (i in 1:nrow(cross_species_cluster_subset)) {
    
    shared_sp <- cross_species_cluster_subset[i, "species"]
    
    gene_type <- c("intact", "pseudo", "both")
    pangenome_type <- c("ultra.cloud", "other.cloud", "soft.core", "shell")
    
    for (g in gene_type) {
      
      for (p in pangenome_type) {
        
        if (cross_species_cluster %in% cluster_pangenome_categories[[shared_sp]][[g]][[p]]) {
          
          category_vec <- c(category_vec, paste(g, p, sep = "|"))
          break
          
        }
      }
    }
  }
  
  names(category_vec) <- cross_species_cluster_subset$species
  
  cross_species_cluster_summary[[cross_species_cluster]] <- category_vec
  
}

# Save RDS containing cross-species cluster information.
saveRDS(object = cross_species_cluster_summary,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cross_species_cluster_summary.rds")



