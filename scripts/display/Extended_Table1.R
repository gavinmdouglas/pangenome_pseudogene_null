rm(list = ls(all.names = TRUE))

### Create table summarizing this information, per species:
###  - Number of genomes
###  - Number of intact genes (within passing clusters) per genome
###  - Number of intact clusters
###  - Number of pseudogenes (within passing clusters) per genome
###  - Number of pseudogene clusters
###  - Mean number of intact clusters per genome
###  - Mean number of pseudogene clusters per genome
###  (Get standard deviations for all mean computations too - just in case they are useful to report!)

# Initialize table.
species <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/focal_species.txt",
                      header = FALSE, stringsAsFactors = FALSE)$V1
summary_means <- data.frame(matrix(NA, nrow = length(species), ncol = 17))
rownames(summary_means) <- species
colnames(summary_means) <- c("Num_genomes",
                             "Total_intact_genes",
                             "Total_pseudo_genes",
                             "Total_intact_clusters",
                             "Total_pseudo_clusters",
                             "Mean_intact_genes",
                             "Mean_intact_clusters",
                             "Mean_pseudo_genes",
                             "Mean_pseudo_clusters",
                             "Mean_intact_percent_cov",
                             "Mean_pseudo_percent_cov",
                             "Sd_intact_genes",
                             "Sd_intact_clusters",
                             "Sd_pseudo_genes",
                             "Sd_pseudo_clusters",
                             "Sd_intact_percent_cov",
                             "Sd_pseudo_percent_cov")

# Read in cluster membership table + definitions of cluster types.
cluster_breakdown <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

cluster_types <- readRDS(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_types.rds")

per_genome_coverages <- readRDS(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/per_genome_element.type_percent_coverages.rds")

# First remove all clusters that are a mix of intact and pseudogene elements.
cluster_breakdown <- cluster_breakdown[-which(cluster_breakdown$cluster %in% cluster_types$both), ]

for (sp in species) {
  cluster_breakdown_sp <- cluster_breakdown[which(cluster_breakdown$species == sp), ]
  cluster_breakdown_sp_intact <- cluster_breakdown_sp[which(cluster_breakdown_sp$cluster %in% cluster_types$intact), ]
  cluster_breakdown_sp_pseudo <- cluster_breakdown_sp[which(cluster_breakdown_sp$cluster %in% cluster_types$intergenic.pseudogene), ]
  
  summary_means[sp, "Num_genomes"] <- length(unique(cluster_breakdown_sp$accession))
  summary_means[sp, "Total_intact_genes"] <- nrow(cluster_breakdown_sp_intact)
  summary_means[sp, "Total_pseudo_genes"] <- nrow(cluster_breakdown_sp_pseudo)
  summary_means[sp, "Total_intact_clusters"] <- length(unique(cluster_breakdown_sp_intact$cluster))
  summary_means[sp, "Total_pseudo_clusters"] <- length(unique(cluster_breakdown_sp_pseudo$cluster))
  summary_means[sp, "Mean_intact_genes"] <- mean(table(cluster_breakdown_sp_intact$accession))
  summary_means[sp, "Mean_pseudo_genes"] <- mean(table(cluster_breakdown_sp_pseudo$accession))
  summary_means[sp, "Mean_intact_percent_cov"] <- mean(per_genome_coverages$intact[[sp]]$genome_percent)
  summary_means[sp, "Mean_pseudo_percent_cov"] <- mean(per_genome_coverages$intergenic.pseudogene[[sp]]$genome_percent)
  
  summary_means[sp, "Sd_intact_genes"] <- sd(table(cluster_breakdown_sp_intact$accession))
  summary_means[sp, "Sd_pseudo_genes"] <- sd(table(cluster_breakdown_sp_pseudo$accession))
  summary_means[sp, "Sd_intact_percent_cov"] <- sd(per_genome_coverages$intact[[sp]]$genome_percent)
  summary_means[sp, "Sd_pseudo_percent_cov"] <- sd(per_genome_coverages$intergenic.pseudogene[[sp]]$genome_percent)
  
  # Additional steps needed to get mean number of clusters.
  cluster_breakdown_sp_intact_clusters <- cluster_breakdown_sp_intact[, c("cluster", "accession")]
  rownames(cluster_breakdown_sp_intact_clusters) <- NULL
  cluster_breakdown_sp_intact_clusters <- cluster_breakdown_sp_intact_clusters[-which(duplicated(cluster_breakdown_sp_intact_clusters)), ]
  summary_means[sp, "Mean_intact_clusters"] <- mean(table(cluster_breakdown_sp_intact_clusters$accession))
  summary_means[sp, "Sd_intact_clusters"] <- sd(table(cluster_breakdown_sp_intact_clusters$accession))
  
  cluster_breakdown_sp_pseudo_clusters <- cluster_breakdown_sp_pseudo[, c("cluster", "accession")]
  rownames(cluster_breakdown_sp_pseudo_clusters) <- NULL
  cluster_breakdown_sp_pseudo_clusters <- cluster_breakdown_sp_pseudo_clusters[-which(duplicated(cluster_breakdown_sp_pseudo_clusters)), ]
  summary_means[sp, "Mean_pseudo_clusters"] <- mean(table(cluster_breakdown_sp_pseudo_clusters$accession))
  summary_means[sp, "Sd_pseudo_clusters"] <- sd(table(cluster_breakdown_sp_pseudo_clusters$accession))
  
}

write.table(x  = summary_means,
            file = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_table_indepth_per_species_breakdown.tsv",
            quote = FALSE,
            row.names = TRUE,
            col.names = NA,
            sep = "\t")

# Note: this table was then formatted in Word and saved manually as a PDF and then converted to JPEG.
