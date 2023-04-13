rm(list = ls(all.names = TRUE))

# Parse CARD annotation and get mapping to clusters.

# Easy for the annotations performed on the both/intact representative clusters. These can simply me re-mapped to the corresponding cluster id.
both_and_intact_representative_seq_CARD <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/card_annot/output/both_and_intact.txt.gz",
                                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

rep2cluster <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/cdhit_workflow_intact_and_intergenic/cluster2rep.tsv.gz",
                          header = FALSE, sep = "\t", row.names = 2, stringsAsFactors = FALSE)

both_and_intact_cluster_CARD <- both_and_intact_representative_seq_CARD[, c("Drug.Class", "Resistance.Mechanism")]
both_and_intact_cluster_CARD$Drug.Class <- gsub("; ", ";", both_and_intact_cluster_CARD$Drug.Class)
both_and_intact_cluster_CARD$Resistance.Mechanism <- gsub("; ", ";", both_and_intact_cluster_CARD$Resistance.Mechanism)
both_and_intact_cluster_CARD$cluster <- rep2cluster[both_and_intact_representative_seq_CARD$ORF_ID, 1]
both_and_intact_cluster_CARD <- both_and_intact_cluster_CARD[, c("cluster", "Drug.Class", "Resistance.Mechanism")]


# Simplifying UniRef-based annotation is trickier, as there are multiple UniRef hits per pseudogene sequence,
# and then those pseudogenes themselves are in clusters.
# Used same approach for regrouping as done for the COG annotations for these data.
pseudogene_uniref_hits_CARD <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/card_annot/output/UniRef_filt_hits.txt.gz",
                                          header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "")

# Simplify to save RAM
pseudogene_uniref_hits_CARD <- pseudogene_uniref_hits_CARD[, c("ORF_ID", "Drug.Class", "Resistance.Mechanism")]

pseudo2uniref <- data.frame(data.table::fread("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/intergenic_uniref_matches.tsv.gz"))

# Subset to pseudo sequences that have at least one annotated uniref hit.
annot_rows_i <- which(pseudo2uniref$uniref_hit %in% pseudogene_uniref_hits_CARD$ORF_ID)
pseudos_w_hit <- unique(pseudo2uniref[annot_rows_i, "intergenic_pseudo_id"])
pseudo2uniref_subset <- pseudo2uniref[which(pseudo2uniref$intergenic_pseudo_id %in% pseudos_w_hit), ]

pseudo2uniref_subset_list <- collapse::rsplit(x = pseudo2uniref_subset$uniref_hit, fl = pseudo2uniref_subset$intergenic_pseudo_id)

pseudo_annot_raw <- parallel::mclapply(X = pseudo2uniref_subset_list,
                                       
                                       FUN = function(uniref_hits) {
                                           
                                         out_df <- data.frame(majority_Drug.Class = "", any_Drug.Class = "",
                                                              majority_Resistance.Mechanism = "", any_Resistance.Mechanism = "")
                                         
                                         hits_annot <- pseudogene_uniref_hits_CARD[pseudogene_uniref_hits_CARD$ORF_ID %in% uniref_hits, ]
                                         
                                         # Return NA only if no hits annotated.
                                         if (nrow(hits_annot) == 0) { return(out_df) }
                                         
                                         hit_annot_Drug.Class <- do.call(c, lapply(hits_annot$Drug.Class, function(x) { strsplit(x = x, split = ";")[[1]] }))
                                         hit_annot_Resistance.Mechanism <- do.call(c, lapply(hits_annot$Resistance.Mechanism, function(x) { strsplit(x = x, split = ";")[[1]] }))
                                        
                                         hit_annot_Drug.Class <- gsub("^ ", "", hit_annot_Drug.Class)
                                         hit_annot_Drug.Class <- gsub(" $", "", hit_annot_Drug.Class)
                                         
                                         hit_annot_Resistance.Mechanism <- gsub("^ ", "", hit_annot_Resistance.Mechanism)
                                         hit_annot_Resistance.Mechanism <- gsub(" $", "", hit_annot_Resistance.Mechanism)
                                         
                                          # Annotate based on majority-rule if possible.
                                         majority_count_required <- ceiling((length(uniref_hits) / 2) + 0.5)
                                         
                                         if (nrow(hits_annot) >= majority_count_required) {
                                           
                                           hit_annot_Drug.Class_tallies <- table(hit_annot_Drug.Class)
                                           majority_Drug.Class_observed_i <- which(hit_annot_Drug.Class_tallies >= majority_count_required)
                                           if (length(majority_Drug.Class_observed_i) > 0) {
                                             out_df$majority_Drug.Class <- paste(unique(names(hit_annot_Drug.Class_tallies)[majority_Drug.Class_observed_i]), collapse = ";")
                                           }
                                           
                                           hit_annot_Resistance.Mechanism_tallies <- table(hit_annot_Resistance.Mechanism)
                                           majority_Resistance.Mechanism_observed_i <- which(hit_annot_Resistance.Mechanism_tallies >= majority_count_required)
                                           if (length(majority_Resistance.Mechanism_observed_i) > 0) {
                                             out_df$majority_Resistance.Mechanism <- paste(unique(names(hit_annot_Resistance.Mechanism_tallies)[majority_Resistance.Mechanism_observed_i]), collapse = ";")
                                           }
                                           
                                         }
                                         
                                         # Also annotate based on "any" annotation approach.
                                         out_df$any_Drug.Class <- paste(sort(unique(hit_annot_Drug.Class)), collapse = ";")
                                         out_df$any_Resistance.Mechanism <- paste(sort(unique(hit_annot_Resistance.Mechanism)), collapse = ";")

                                         return(out_df)
                                       },
                                       mc.cores = 14)

pseudo_individual_annot <- do.call(rbind, pseudo_annot_raw)


# Then do similar regrouping step, but to get annotations at the cluster level.
cluster_types <- readRDS(file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_types.rds")

# Read in cluster membership table + definitions of cluster types.
cluster_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/summary/cluster_member_breakdown.tsv.gz",
                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

# Sanity check that all pseudo ids are in the cluster breakdown table:
length(which(! rownames(pseudo_individual_annot) %in% rownames(cluster_breakdown)))

# The above was non-zero, so removed rows missing in cluster breakdown table
# (which can happen if they were only present in genomes that were filtered out).
pseudo_individual_annot <- pseudo_individual_annot[which(rownames(pseudo_individual_annot) %in% rownames(cluster_breakdown)), ]

# Regroup intergenic elements to individual clusters.
# Note that in this case the "any" functions are only reported if >50% of underlying pseudogenes
# have that function marked as "any" based on the preceding UniRef hit analysis.
pseudo_individual_annot$cluster <- cluster_breakdown[rownames(pseudo_individual_annot), "cluster"]

cluster2pseudo_list <- collapse::rsplit(x = rownames(pseudo_individual_annot), fl = pseudo_individual_annot$cluster)

pseudo_cluster_annot_raw <- parallel::mclapply(X = names(cluster2pseudo_list),
                                           
                                           FUN = function(cluster_id) {
                                             
                                             out_df <- data.frame(majority_Drug.Class = "", any_Drug.Class = "",
                                                                  majority_Resistance.Mechanism = "", any_Resistance.Mechanism = "")
                                             
                                             hits_annot <- pseudo_individual_annot[cluster2pseudo_list[[cluster_id]], ]
                                             
                                             # Return error if number of hits returned does not match expected.
                                             if (nrow(hits_annot) != length(cluster2pseudo_list[[cluster_id]])) { stop("Mismatch with number of observed and expected hits") }
                                             
                                             # First process Drug.Class data.
                                             hit_annot_majority_Drug.Class <- do.call(c, lapply(hits_annot$majority_Drug.Class, function(x) { strsplit(x = x, split = ";")[[1]] }))
                                             hit_annot_any_Drug.Class <- do.call(c, lapply(hits_annot$any_Drug.Class, function(x) { strsplit(x = x, split = ";")[[1]] }))
                                             
                                             hit_annot_majority.Drug.Class_tallies <- table(hit_annot_majority_Drug.Class)
                                             hit_annot_any.Drug.Class_tallies <- table(hit_annot_any_Drug.Class)
                                             
                                             # Annotate based on majority-rule if possible.
                                             majority_count_required <- ceiling((nrow(hits_annot) / 2) + 0.5)
                                             
                                             hit_annot_majority.Drug.Class_tallies_observed_i <- which(hit_annot_majority.Drug.Class_tallies >= majority_count_required)
                                             hit_annot_any.Drug.Class_tallies_observed_i <- which(hit_annot_any.Drug.Class_tallies >= majority_count_required)
                                             
                                             if (length(hit_annot_majority.Drug.Class_tallies_observed_i) > 0) {
                                               out_df$majority_Drug.Class <- paste(names(hit_annot_majority.Drug.Class_tallies)[hit_annot_majority.Drug.Class_tallies_observed_i], collapse = ";")
                                             }

                                             if (length(hit_annot_any.Drug.Class_tallies_observed_i) > 0) {
                                               out_df$any_Drug.Class <- paste(names(hit_annot_any.Drug.Class_tallies)[hit_annot_any.Drug.Class_tallies_observed_i], collapse = ";")
                                             }
                                            
                                             
                                             # First process Resistance.Mechanism data.
                                             hit_annot_majority_Resistance.Mechanism <- do.call(c, lapply(hits_annot$majority_Resistance.Mechanism, function(x) { strsplit(x = x, split = ";")[[1]] }))
                                             hit_annot_any_Resistance.Mechanism <- do.call(c, lapply(hits_annot$any_Resistance.Mechanism, function(x) { strsplit(x = x, split = ";")[[1]] }))
                                             
                                             hit_annot_majority.Resistance.Mechanism_tallies <- table(hit_annot_majority_Resistance.Mechanism)
                                             hit_annot_any.Resistance.Mechanism_tallies <- table(hit_annot_any_Resistance.Mechanism)
                                             
                                             # Annotate based on majority-rule if possible.
                                             majority_count_required <- ceiling((nrow(hits_annot) / 2) + 0.5)
                                             
                                             hit_annot_majority.Resistance.Mechanism_tallies_observed_i <- which(hit_annot_majority.Resistance.Mechanism_tallies >= majority_count_required)
                                             hit_annot_any.Resistance.Mechanism_tallies_observed_i <- which(hit_annot_any.Resistance.Mechanism_tallies >= majority_count_required)
                                             
                                             if (length(hit_annot_majority.Resistance.Mechanism_tallies_observed_i) > 0) {
                                               out_df$majority_Resistance.Mechanism <- paste(names(hit_annot_majority.Resistance.Mechanism_tallies)[hit_annot_majority.Resistance.Mechanism_tallies_observed_i], collapse = ";")
                                             }
                                             
                                             if (length(hit_annot_any.Resistance.Mechanism_tallies_observed_i) > 0) {
                                               out_df$any_Resistance.Mechanism <- paste(names(hit_annot_any.Resistance.Mechanism_tallies)[hit_annot_any.Resistance.Mechanism_tallies_observed_i], collapse = ";")
                                             }
                                              
                                             return(out_df)

                                           },
                                           mc.cores = 14)

pseudo_cluster_annot <- do.call(rbind, pseudo_cluster_annot_raw)
pseudo_cluster_annot$cluster <- names(cluster2pseudo_list)
pseudo_cluster_annot_both <- pseudo_cluster_annot[which(pseudo_cluster_annot$cluster %in% cluster_types$both), ]
pseudo_cluster_annot_pseudo.only <- pseudo_cluster_annot[which(pseudo_cluster_annot$cluster %in% cluster_types$intergenic.pseudogene), ]

pseudo_cluster_annot_both <- pseudo_cluster_annot_both[, c("cluster", "majority_Drug.Class", "any_Drug.Class",
                                                           "majority_Resistance.Mechanism", "any_Resistance.Mechanism")]

pseudo_cluster_annot_pseudo.only <- pseudo_cluster_annot_pseudo.only[, c("cluster", "majority_Drug.Class", "any_Drug.Class",
                                                                         "majority_Resistance.Mechanism", "any_Resistance.Mechanism")]

pseudo_individual_annot <- pseudo_individual_annot[, c("cluster", "majority_Drug.Class", "any_Drug.Class",
                                                       "majority_Resistance.Mechanism", "any_Resistance.Mechanism")]

# Save all output annotations at cluster level in RDS.
saveRDS(object = list(both_and_intact_cluster_CARD = both_and_intact_cluster_CARD,
                       pseudo_individual_element_annot = pseudo_individual_annot,
                       pseudo_cluster_CARD = pseudo_cluster_annot_pseudo.only,
                       both_uniref.based_cluster_CARD = pseudo_cluster_annot_both),
        file = "/data1/gdouglas/projects/accessory_vs_pseudogene/summary/CARD_cluster_annot_by_type.rds")
     
