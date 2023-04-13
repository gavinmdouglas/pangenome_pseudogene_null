# Parse EGG-NOG output for UniRef ref matches that hit intergenic pseudogenes (based on pseudofinder).

rm(list = ls(all.names = TRUE))

source("/home/gdouglas/scripts/accessory_vs_pseudogenes/R_scripts/functions.R")

COG2catgory_uniq <- readRDS("/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")

ref_eggNOG_annot <- data.table::fread("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/eggnog_uniref/eggnog_mapper_annot/uniref_annot.emapper.annotations.tsv.gz")

ref_eggNOG_annot$all_COG <- sapply(ref_eggNOG_annot$eggNOG_OGs, function(x) { paste(grep("^COG", unique(gsub("@.*$", "", strsplit(x, ",")[[1]])), value = TRUE), collapse = ",") })

# Remove hits with no COG annotation.
ref_eggNOG_annot <- ref_eggNOG_annot[which(ref_eggNOG_annot$all_COG != ""), ]

# Simplify to save RAM
ref_eggNOG_annot <- ref_eggNOG_annot[, c("query", "all_COG")]

pseudo2uniref <- data.table::fread("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_uniref_matches.tsv.gz")

pseudo2uniref_list <- collapse::rsplit(x = pseudo2uniref$uniref_hit, fl = pseudo2uniref$intergenic_pseudo_id)

pseudo_annot_raw <- parallel::mclapply(X = pseudo2uniref_list,

                                       FUN = function(uniref_hits) {
                                         
                                         out_df <- data.frame(majority_COG = "", majority_COG_category = "", any_COG = "", any_COG_category = "")
                                         
                                         hits_annot <- ref_eggNOG_annot[ref_eggNOG_annot$query %in% uniref_hits, ]
                                         
                                         # Return NA only if no hits annotated.
                                         if (nrow(hits_annot) == 0) { return(out_df) }
                                         
                                         hit_annot_COG <- do.call(c, lapply(hits_annot$all_COG, function(x) { strsplit(x = x, split = ",")[[1]] }))
                                         
                                         # Annotate based on majority-rule if possible.
                                         majority_count_required <- ceiling((length(uniref_hits) / 2) + 0.5)

                                         if (nrow(hits_annot) >= majority_count_required) {
                                         
                                           hit_annot_COG_tallies <- table(hit_annot_COG)
                                           
                                           majority_observed_i <- which(hit_annot_COG_tallies >= majority_count_required)
                                           
                                           if (length(majority_observed_i) > 0) {
                                             out_df$majority_COG <- paste(names(hit_annot_COG_tallies)[majority_observed_i], collapse = ",")
                                             raw_majority_COG_categories <- COG2catgory_uniq[names(hit_annot_COG_tallies)[majority_observed_i], "category"]
                                             out_df$majority_COG_category <- paste(unique(sort(do.call(c, lapply(raw_majority_COG_categories, function(x) { strsplit(x = x, split = ",")[[1]] })))), collapse = ",")
                                           }
                                          
                                         }
                                         
                                         
                                         # Also annotate based on "any" annotation approach.
                                         unique_COG_annot <- sort(unique(hit_annot_COG))
                                         out_df$any_COG <- paste(unique_COG_annot, collapse = ",")
                                         out_df$any_COG_category <- paste(unique(sort(do.call(c, lapply(COG2catgory_uniq[unique_COG_annot, "category"], function(x) { strsplit(x = x, split = ",")[[1]] })))), collapse = ",")

                                         return(out_df)
                                       },
                                       mc.cores = 14)
  


pseudo_annot <- do.call(rbind, pseudo_annot_raw)

saveRDS(object = pseudo_annot,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intergenic_pseudo_uniref_based_COG_annot.rds")
