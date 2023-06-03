rm(list = ls(all.names = TRUE))

# Summarize AIC, both absolute and normalized, and variance explained by each random effect variable for final models.

library(ggplot2)

glmm_fit.info_RAW <- list()
glmm_final_varcor_RAW <- list()

# First, compare the model fits based on varying parameter complexity.
for (partition in c("shell", "other.cloud", "ultra.cloud")) {
  
  glmm_summaries <- readRDS(paste("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/glmm_output/", partition, "_model_output_summaries.rds", sep = ""))
  glmm_summaries_additional <- readRDS(paste("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/glmm_output/", partition, "_model_output_summaries_additional.rds", sep = ""))
  glmm_summaries[["species.only"]] <- glmm_summaries_additional$species_only
  
  
  random_effect_sd <- sapply(glmm_summaries$COG.only_species.interaction_redundant.both.interaction$varcor$cond,
                             function(x) { attr(x, which = "stddev") })
  
  glmm_final_varcor_RAW[[partition]] <- data.frame(partition = partition,
                                                   random_effect = names(glmm_summaries$COG.only_species.interaction_redundant.both.interaction$varcor$cond),
                                                   sd = random_effect_sd) 
  
  glmm_fit.info_RAW[[partition]] <- data.frame(matrix(NA, nrow = length(glmm_summaries), ncol = 5))
  colnames(glmm_fit.info_RAW[[partition]] ) <- names(glmm_summaries$COG.only$AICtab)
  rownames(glmm_fit.info_RAW[[partition]] ) <- names(glmm_summaries)
  
  for (category in names(glmm_summaries)) {
    glmm_fit.info_RAW[[partition]][category, ] <- as.numeric(glmm_summaries[[category]]$AICtab)
  }
  
  glmm_fit.info_RAW[[partition]]$min_by_AIC <- glmm_fit.info_RAW[[partition]]$AIC - min(glmm_fit.info_RAW[[partition]]$AIC)
  glmm_fit.info_RAW[[partition]]$min_by_AIC_norm <- glmm_fit.info_RAW[[partition]]$min_by_AIC / max(glmm_fit.info_RAW[[partition]]$min_by_AIC)
  
  glmm_fit.info_RAW[[partition]]$partition <- partition
  glmm_fit.info_RAW[[partition]]$model <- rownames(glmm_fit.info_RAW[[partition]])

}

glmm_fit.info <- do.call(rbind, glmm_fit.info_RAW)
rownames(glmm_fit.info) <- NULL

glmm_fit.info$partition_clean <- NA
glmm_fit.info$partition_clean[which(glmm_fit.info$partition == "shell")] <- "Shell"
glmm_fit.info$partition_clean[which(glmm_fit.info$partition == "other.cloud")] <- "Other-rare"
glmm_fit.info$partition_clean[which(glmm_fit.info$partition == "ultra.cloud")] <- "Ultra-rare"
glmm_fit.info$partition_clean <- factor(glmm_fit.info$partition_clean, levels = c("Ultra-rare", "Other-rare", "Shell"))

glmm_fit.info$model <- gsub("COG.only_", "COG_", glmm_fit.info$model)
glmm_fit.info$model <- factor(glmm_fit.info$model, levels = rev(c("species.only",
                                                                  "COG.only",
                                                                  "COG_species",
                                                                  "COG_species.interaction",
                                                                  "COG_species.interaction_redundant",
                                                                  "COG_species.interaction_redundant.species.interaction",
                                                                  "COG_species.interaction_redundant.COG.interaction",
                                                                  "COG_species.interaction_redundant.both.interaction")))

glmm_fit.info$formula <- NA
glmm_fit.info[which(glmm_fit.info$model == "species.only"), "formula"] <- "pseudogene ~ (1 | species)"
glmm_fit.info[which(glmm_fit.info$model == "COG.only"), "formula"] <- "pseudogene ~ COG-category"
glmm_fit.info[which(glmm_fit.info$model == "COG_species"), "formula"] <- "pseudogene ~ COG-category + (1 | species)"
glmm_fit.info[which(glmm_fit.info$model == "COG_species.interaction"), "formula"] <- "pseudogene ~ COG-category + (1 | species) +\n(1 | COG_category : species)"
glmm_fit.info[which(glmm_fit.info$model == "COG_species.interaction_redundant"), "formula"] <- "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG_category : species)"
glmm_fit.info[which(glmm_fit.info$model == "COG_species.interaction_redundant.species.interaction"), "formula"] <- "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG-category : species) + (1 | non-redundant-status : species)"
glmm_fit.info[which(glmm_fit.info$model == "COG_species.interaction_redundant.COG.interaction"), "formula"] <- "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG-category : species) +\n(1 | non-redundant-status : COG-category)"
glmm_fit.info[which(glmm_fit.info$model == "COG_species.interaction_redundant.both.interaction"), "formula"] <- "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG-category : species) +\n(1 | non-redundant-status : species) +\n(1 | non-redundant-status : COG-category)"

glmm_fit.info$formula <- factor(glmm_fit.info$formula, levels = rev(c("pseudogene ~ (1 | species)",
                                                                      "pseudogene ~ COG-category",
                                                                      "pseudogene ~ COG-category + (1 | species)",
                                                                      "pseudogene ~ COG-category + (1 | species) +\n(1 | COG_category : species)",
                                                                      "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG_category : species)",
                                                                      "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG-category : species) + (1 | non-redundant-status : species)",
                                                                      "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG-category : species) +\n(1 | non-redundant-status : COG-category)",
                                                                      "pseudogene ~ COG-category + non-redundant-status + (1 | species) +\n(1 | COG-category : species) +\n(1 | non-redundant-status : species) +\n(1 | non-redundant-status : COG-category)")))

# Add raw AIC values to this barplot instead of a separate table.
glmm_fit.info$AIC_round <- format(round(glmm_fit.info$AIC, digits=2), nsmall = 1)

glmm_fit.info$hjust_setting <- -0.05
glmm_fit.info$text_colour <- "black"
glmm_fit.info[which(glmm_fit.info$min_by_AIC_norm >= 0.8), "hjust_setting"] <- 1.075
glmm_fit.info[which(glmm_fit.info$min_by_AIC_norm >= 0.8), "text_colour"] <- "white"
norm_AIC_barplot <- ggplot(data = glmm_fit.info, aes(x = min_by_AIC_norm, y = formula)) +
                           geom_bar(stat = "identity") +
                           geom_text(aes(label=AIC_round), colour = glmm_fit.info$text_colour, hjust = glmm_fit.info$hjust_setting) +
                           facet_wrap(partition_clean ~ ., scales = "free_x") +
                           xlab("Normalized AIC: (AIC - min[AIC]) / max(AIC - min(AIC))") +
                           ylab("") +
                           theme_bw()

ggsave(plot = norm_AIC_barplot,
       filename = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_indepth_GLMM_AIC_comparison.pdf",
       device = "pdf", width = 13, height = 10, units = "in", dpi = 400)


## Summarize variance explained by each random effect.
## Optional plot that could also be included.
# glmm_final_varcor <- do.call(rbind, glmm_final_varcor_RAW)
# rownames(glmm_final_varcor) <- NULL
# 
# glmm_final_varcor$random_effect_formula <- NA
# glmm_final_varcor$random_effect_formula[which(glmm_final_varcor$random_effect == "species")] <- "(1 | species)"
# glmm_final_varcor$random_effect_formula[which(glmm_final_varcor$random_effect == "COG_category:species")] <- "(1 | COG-category : species)"
# glmm_final_varcor$random_effect_formula[which(glmm_final_varcor$random_effect == "NO_redundant_intact_COG:species")] <- "(1 | non-redundant-status : species)"
# 
# glmm_final_varcor$random_effect_formula <- factor(glmm_final_varcor$random_effect_formula,
#                                            levels = rev(c("(1 | species)", "(1 | COG-category : species)", "(1 | non-redundant-status : species)")))
# 
# glmm_final_varcor$partition_clean <- NA
# glmm_final_varcor$partition_clean[which(glmm_final_varcor$partition == "shell")] <- "Shell"
# glmm_final_varcor$partition_clean[which(glmm_final_varcor$partition == "other.cloud")] <- "Other-rare"
# glmm_final_varcor$partition_clean[which(glmm_final_varcor$partition == "ultra.cloud")] <- "Ultra-rare"
# glmm_final_varcor$partition_clean <- factor(glmm_final_varcor$partition_clean, levels = c("Ultra-rare", "Other-rare", "Shell"))
# 
# ggplot(data = glmm_final_varcor, aes(x = sd, y = random_effect_formula)) +
#   geom_bar(stat = "identity") +
#   facet_wrap(partition_clean ~ ., scales = "free_x") +
#   xlab("Standard deviation") +
#   ylab("Random effect") +
#   theme_bw()
# 
