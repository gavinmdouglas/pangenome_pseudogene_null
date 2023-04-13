rm(list = ls(all.names = TRUE))

COG2catgory_uniq <- readRDS("/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")

# Read in cluster membership table + definitions of cluster types.
cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

cluster_types <- readRDS(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_types.rds")


# Read in COG annotations for intergenic elements based on UniRef hit annotations.
intergenic_uniref_COG <- readRDS("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intergenic_pseudo_uniref_based_COG_annot.rds")

# Remove intergenic elements that were filtered out (and so are not in cluster table).
intergenic_uniref_COG <- intergenic_uniref_COG[intersect(rownames(intergenic_uniref_COG), rownames(cluster_breakdown)), ]

# Regroup intergenic elements to individual clusters.
intergenic_uniref_COG$cluster <- cluster_breakdown[rownames(intergenic_uniref_COG), "cluster"]

clusters_w_majority.COG <- grep("^COG", intergenic_uniref_COG$majority_COG)

cluster2pseudo_list <- collapse::rsplit(x = rownames(intergenic_uniref_COG), fl = intergenic_uniref_COG$cluster)

intergenic_annot_raw <- parallel::mclapply(X = names(cluster2pseudo_list),

                                           FUN = function(cluster_id) {

                                           out_df <- data.frame(cluster = cluster_id, majority_COG = "", majority_COG_category = "", any_COG = "", any_COG_category = "")

                                           hits_annot <- intergenic_uniref_COG[cluster2pseudo_list[[cluster_id]], ]

                                           # Return error if number of hits returned does not match expected.
                                           if (nrow(hits_annot) != length(cluster2pseudo_list[[cluster_id]])) { stop("Mismatch with number of observed and expected hits") }

                                           hit_annot_majority.COG <- do.call(c, lapply(hits_annot$majority_COG, function(x) { strsplit(x = x, split = ",")[[1]] }))
                                           hit_annot_any_COG <- do.call(c, lapply(hits_annot$any_COG, function(x) { strsplit(x = x, split = ",")[[1]] }))

                                           hit_annot_majority.COG_tallies <- table(hit_annot_majority.COG)
                                           hit_annot_any.COG_tallies <- table(hit_annot_any_COG)

                                           majority_count_required <- ceiling((nrow(hits_annot) / 2) + 0.5)

                                           hit_annot_majority.COG_tallies_majority_observed_i <- which(hit_annot_majority.COG_tallies >= majority_count_required)
                                           hit_annot_any.COG_tallies_majority_observed_i <- which(hit_annot_any.COG_tallies >= majority_count_required)

                                           if (length(hit_annot_majority.COG_tallies_majority_observed_i) > 0) {
                                             out_df$majority_COG <- paste(names(hit_annot_majority.COG_tallies)[hit_annot_majority.COG_tallies_majority_observed_i], collapse = ",")
                                             raw_majority_COG_categories <- COG2catgory_uniq[names(hit_annot_majority.COG_tallies)[hit_annot_majority.COG_tallies_majority_observed_i], "category"]
                                             out_df$majority_COG_category <- paste(unique(sort(do.call(c, lapply(raw_majority_COG_categories, function(x) { strsplit(x = x, split = ",")[[1]] })))), collapse = ",")
                                           }

                                           if (length(hit_annot_any.COG_tallies_majority_observed_i) > 0) {
                                             out_df$any_COG <- paste(names(hit_annot_any.COG_tallies)[hit_annot_any.COG_tallies_majority_observed_i], collapse = ",")
                                             raw_any_COG_categories <- COG2catgory_uniq[names(hit_annot_any.COG_tallies)[hit_annot_any.COG_tallies_majority_observed_i], "category"]
                                             out_df$any_COG_category <- paste(unique(sort(do.call(c, lapply(raw_any_COG_categories, function(x) { strsplit(x = x, split = ",")[[1]] })))), collapse = ",")
                                           }
                                           
                                           return(out_df)

                                         },
                                         mc.cores = 14)

intergenic_annot <- do.call(rbind, intergenic_annot_raw)

intergenic_annot_both <- intergenic_annot[which(intergenic_annot$cluster %in% cluster_types$both), ]
intergenic_annot_pseudo.only <- intergenic_annot[which(intergenic_annot$cluster %in% cluster_types$intergenic.pseudogene), ]

all_intergenic_annot <- list(pseudo = intergenic_annot_pseudo.only,
                             both_intact_and_pseudo_uniref_approach = intergenic_annot_both)

saveRDS(object = all_intergenic_annot,
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intergenic_pseudo_cluster_COG_annot.rds")
