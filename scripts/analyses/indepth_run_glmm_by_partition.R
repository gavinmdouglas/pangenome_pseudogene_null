rm(list = ls(all.names = TRUE))

library(glmmTMB)

num_cores <- 5

element_info <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/element_glmm_input.tsv.gz",
                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Remove rows corresponding to really rare COG categories
to_ignore <- c("A", "B", "Y", "Z")
to_ignore_category_row_i <- which(element_info$COG_category %in% to_ignore)
#length(to_ignore_category_row_i)
element_info <- element_info[-to_ignore_category_row_i, ]

element_info$pseudogenized <- FALSE
element_info[grep("pseudo", rownames(element_info)), "pseudogenized"] <- TRUE

# Set K ("Transcription") as reference COG category.
COG_category_level_order <- sort(unique(element_info$COG_category))
COG_category_level_order <- c("K", COG_category_level_order[-which(COG_category_level_order == "K")])
element_info$COG_category <- factor(element_info$COG_category, levels = COG_category_level_order)

# Reverse redundant_intact_COG, as the majority of elements have a redundant COG, so this makes more sense as the reference.
element_info$NO_redundant_intact_COG <- ! element_info$redundant_intact_COG

### Code used to figure out what COG category reference would make most sense (based on abundance).
# breakdown <- list()
# for (partition in c("shell", "other.cloud", "ultra.cloud")) {
#   element_info_subset <- element_info[which(element_info$partition == partition), ]
#   breakdown[[partition]] <- table(element_info_subset$COG_category)
# }
# sapply(breakdown, sort, decreasing = TRUE)
### Based on this "K" is the most consistently abundant COG category across all three partitions
### (#3 in other and shell, and #4 in ultra). 

for (partition in c("shell", "other.cloud", "ultra.cloud")) {
  
  element_info_subset <- element_info[which(element_info$partition == partition), ]
  
  model_output <- list()
  
  model_output[["COG.only"]] <- glmmTMB(formula = pseudogenized ~ COG_category,
                                        data = element_info_subset,
                                        family = "binomial",
                                        control = glmmTMBControl(optimizer = nlminb,
                                                                 parallel = num_cores,
                                                                 profile = TRUE,
                                                                 optCtrl = list(iter.max = 1000,
                                                                                eval.max = 1000)))
  model_output[["COG.only"]]$frame <- NULL
  
  model_output[["COG.only_species"]] <- glmmTMB(formula = pseudogenized ~ COG_category + (1 | species),
                                           data = element_info_subset,
                                           family = "binomial",
                                           control = glmmTMBControl(optimizer = nlminb,
                                                                    parallel = num_cores,
                                                                    profile = TRUE,
                                                                    optCtrl = list(iter.max = 1000,
                                                                                   eval.max = 1000)))
  model_output[["COG.only_species"]]$frame <- NULL
  
  model_output[["COG.only_species.interaction"]] <- glmmTMB(formula = pseudogenized ~ COG_category + (1 | species) + (1 | COG_category:species),
                                                       data = element_info_subset,
                                                       family = "binomial",
                                                       control = glmmTMBControl(optimizer = nlminb,
                                                                                parallel = num_cores,
                                                                                profile = TRUE,
                                                                                optCtrl = list(iter.max = 1000,
                                                                                               eval.max = 1000)))
  model_output[["COG.only_species.interaction"]]$frame <- NULL
  
  
  model_output[["COG.only_species.interaction_redundant"]] <- glmmTMB(formula = pseudogenized ~ COG_category + NO_redundant_intact_COG + (1 | species) + (1 | COG_category:species),
                                                                   data = element_info_subset,
                                                                   family = "binomial",
                                                                   control = glmmTMBControl(optimizer = nlminb,
                                                                                            parallel = num_cores,
                                                                                            profile = TRUE,
                                                                                            optCtrl = list(iter.max = 1000,
                                                                                                           eval.max = 1000)))
  model_output[["COG.only_species.interaction_redundant"]]$frame <- NULL
  
  
  model_output[["COG.only_species.interaction_redundant.COG.interaction"]] <- glmmTMB(formula = pseudogenized ~ COG_category + NO_redundant_intact_COG + COG_category:NO_redundant_intact_COG + (1 | species) + (1 | COG_category:species),
                                                                                 data = element_info_subset,
                                                                                 family = "binomial",
                                                                                 control = glmmTMBControl(optimizer = nlminb,
                                                                                                          parallel = num_cores,
                                                                                                          profile = TRUE,
                                                                                                          optCtrl = list(iter.max = 1000,
                                                                                                                         eval.max = 1000)))
  model_output[["COG.only_species.interaction_redundant.COG.interaction"]]$frame <- NULL
  
  
  model_output[["COG.only_species.interaction_redundant.species.interaction"]] <- glmmTMB(formula = pseudogenized ~ COG_category + NO_redundant_intact_COG + (1 | species) + (1 | COG_category:species) + (1 | NO_redundant_intact_COG:species),
                                                                                 data = element_info_subset,
                                                                                 family = "binomial",
                                                                                 control = glmmTMBControl(optimizer = nlminb,
                                                                                                          parallel = num_cores,
                                                                                                          profile = TRUE,
                                                                                                          optCtrl = list(iter.max = 1000,
                                                                                                                         eval.max = 1000)))
  model_output[["COG.only_species.interaction_redundant.species.interaction"]]$frame <- NULL
  
  model_output[["COG.only_species.interaction_redundant.both.interaction"]] <- glmmTMB(formula = pseudogenized ~ COG_category + NO_redundant_intact_COG + COG_category:NO_redundant_intact_COG + (1 | species) + (1 | COG_category:species) + (1 | NO_redundant_intact_COG:species),
                                                                                     data = element_info_subset,
                                                                                     family = "binomial",
                                                                                     control = glmmTMBControl(optimizer = nlminb,
                                                                                                              parallel = num_cores,
                                                                                                              profile = TRUE,
                                                                                                              optCtrl = list(iter.max = 1000,
                                                                                                                             eval.max = 1000)))
  model_output[["COG.only_species.interaction_redundant.both.interaction"]]$frame <- NULL
  
  model_summaries <- list()
  for (model_obj in names(model_output)) {
    model_summaries[[model_obj]] <- summary(model_output[[model_obj]])
  }
  
  saveRDS(object = model_output,
          file = paste("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/", partition, "_model_outputs.rds", sep = ""))
  
  saveRDS(object = model_summaries,
          file = paste("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/", partition, "_model_output_summaries.rds", sep = ""))

}
