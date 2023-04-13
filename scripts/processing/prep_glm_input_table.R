rm(list = ls(all.names = TRUE))

source("/home/gdouglas/scripts/accessory_vs_pseudogenes/R_scripts/functions.R")

cluster_filt_lengths <- read.table(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_filt_lengths_and_additional.tsv.gz",
                                   header = TRUE, sep = "\t", row.names = 1)

cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

# Restrict to clusters that remained after length filters.
cluster_breakdown <- cluster_breakdown[which(cluster_breakdown$cluster %in% rownames(cluster_filt_lengths)), ]

# Read in COG annot information
COG_category_annot <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/all_filt_cluster_COG_annot.rds")

# Read in COG category db info.
COG_category_descrip <- read.table("/data1/gdouglas/db/COG_definitions/COG_category_descrip.tsv",
                                   header = FALSE, sep="\t", stringsAsFactors = FALSE, row.names = 1)

to_ignore <- c("A", "B", "Y", "Z")

cluster_pangenome_categories <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_pangenome_categories.rds")

# Restrict cluster breakdown info to COG-annotated clusters only.
cluster_breakdown <- cluster_breakdown[which(cluster_breakdown$cluster %in% rownames(COG_category_annot)), ]

cluster_breakdown$COG_category <- COG_category_annot[cluster_breakdown$cluster, "COG_category"]

redundant_intact_COG_RAW <- list()
species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt", sep = "\t")$V1
for (sp in species) {
  redundant_intact_COG_RAW[[sp]] <- read.table(paste("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/pseudo_compensation_breakdowns/",
                                                            sp,
                                                            "_redundant_COG_breakdown.tsv.gz",
                                                            sep = ""),
                                                      sep = "\t", stringsAsFactors = FALSE, header = TRUE)
}

redundant_intact_COG <- do.call(rbind, redundant_intact_COG_RAW)

cluster_breakdown$partition <- NA
cluster_breakdown$redundant_intact_COG <- NA

cluster_breakdown[redundant_intact_COG$element, "partition"] <- redundant_intact_COG$partition
cluster_breakdown[redundant_intact_COG$element, "redundant_intact_COG"] <- redundant_intact_COG$redundant_intact_COG


# Many rows were missing from the redundant_intact_COG tables - should be because they are soft.core elements,
# which are not of interest for this analysis.
# Run quick sanity check to confirm this, and remove these rows.
NA_row_i <- which(is.na(cluster_breakdown$redundant_intact_COG) | is.na(cluster_breakdown$partition))
cluster_breakdown_NAs <- cluster_breakdown[NA_row_i, ]
cluster_breakdown <- cluster_breakdown[-NA_row_i, ]

cluster_breakdown_NAs_unique_clusters <- unique(cluster_breakdown_NAs$cluster)

# Get set of all soft.core clusters and accessory gene "both" clusters.
all_soft.core <- as.character()
accessory_both <- as.character()
for (sp in species) {
  for (type in c("pseudo", "intact", "both")) {
    all_soft.core <- c(all_soft.core, cluster_pangenome_categories[[sp]][[type]]$soft.core)
  }
  accessory_both <- c(accessory_both, cluster_pangenome_categories[[sp]][["both"]]$ultra.cloud)
  accessory_both <- c(accessory_both, cluster_pangenome_categories[[sp]][["both"]]$other.cloud)
  accessory_both <- c(accessory_both, cluster_pangenome_categories[[sp]][["both"]]$shell)
  
}

# Indeed all are soft-core, so can ignore them. Note that during the first past, all "both"-classified elements
# were NA as well, but this is no longer the case after a bug fix.
length(which(! cluster_breakdown_NAs_unique_clusters %in% all_soft.core))
length(which(cluster_breakdown_NAs_unique_clusters %in% accessory_both))

# Remove columns not needed for GLM, to save space + time.
cluster_breakdown <- cluster_breakdown[, -which(colnames(cluster_breakdown) %in% c("contig", "cluster", "accession"))]

# Duplicate rows with multiple COG categories, in order to get single category per line.
cluster_breakdown_split <- split_multi_category_rows(in_df = cluster_breakdown,
                                                     category_col = "COG_category",
                                                     delimiter = ",",
                                                     num_cores = 20)

write.table(x = cluster_breakdown_split,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/element_glm_input.tsv",
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = NA)
