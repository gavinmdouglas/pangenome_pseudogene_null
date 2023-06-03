rm(list = ls(all.names = TRUE))

library(ggplot2)

glmm_fit.info_RAW <- list()
glmm_final_summaries_RAW <- list()

# Compare the model fits based on varying parameter complexity.
# Restricted to shell and Other-rare
for (partition in c("shell", "other.cloud")) {

  glmm_summaries <- readRDS(paste("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/", partition, "_model_output_summaries.rds", sep = ""))

  glmm_final_summaries_RAW[[partition]] <- data.frame(glmm_summaries$COG.only_species.interaction_redundant.both.interaction$coefficients$cond)
  glmm_final_summaries_RAW[[partition]]$variable <- rownames(glmm_final_summaries_RAW[[partition]])
  glmm_final_summaries_RAW[[partition]]$partition <- partition

}

glmm_final_summaries <- do.call(rbind, glmm_final_summaries_RAW)
rownames(glmm_final_summaries) <- NULL

glmm_final_summaries$Type <- NA
glmm_final_summaries[which(glmm_final_summaries$variable == "(Intercept)"), "Type"] <- "Intercept"
glmm_final_summaries[which(glmm_final_summaries$variable == "(Intercept)"), "variable"] <- "Intercept"

glmm_final_summaries[grep("NO_redundant", glmm_final_summaries$variable), "Type"] <- "Non-redundant"

COG_category_fixed_effect_i <- intersect(grep("NO_redundant", glmm_final_summaries$variable, invert = TRUE),
                                         grep("COG_category", glmm_final_summaries$variable))
glmm_final_summaries[COG_category_fixed_effect_i, "Type"] <- "COG category"


glmm_final_summaries$variable <- gsub("COG_category", "", glmm_final_summaries$variable)
glmm_final_summaries$variable <- gsub("NO_redundant_intact_COGTRUE", "Non-redundant", glmm_final_summaries$variable)

glmm_final_summaries$Sig <- "No"
glmm_final_summaries$Sig[which(glmm_final_summaries$Pr...z.. < 0.05)] <- "Yes"

# Make combined non-redundant and interaction estimate.
for (partition_level in unique(glmm_final_summaries$partition)) {
  glmm_final_summaries[which(glmm_final_summaries$partition == partition_level & glmm_final_summaries$Type == "Non-redundant"), "Estimate"] <-
    glmm_final_summaries[which(glmm_final_summaries$partition == partition_level & glmm_final_summaries$Type == "Non-redundant"), "Estimate"] +
    glmm_final_summaries[which(glmm_final_summaries$partition == partition_level & glmm_final_summaries$variable == "Non-redundant"), "Estimate"]
}

glmm_final_summaries <- glmm_final_summaries[-which(glmm_final_summaries$variable == "Non-redundant"), ]
glmm_final_summaries$variable <- gsub(":Non-redundant", " (Non-redundant)", glmm_final_summaries$variable)

non_redundant_category_string <- "Non-redundant\n(and interaction)"
glmm_final_summaries$Type[which(glmm_final_summaries$Type == "Non-redundant")] <- non_redundant_category_string

glmm_final_summaries$Type <- factor(glmm_final_summaries$Type, levels = c("Intercept",
                                                                          "COG category",
                                                                          non_redundant_category_string))

glmm_final_summaries$partition_clean <- glmm_final_summaries$partition
glmm_final_summaries$partition_clean[which(glmm_final_summaries$partition == "shell")] <- "Shell"
glmm_final_summaries$partition_clean[which(glmm_final_summaries$partition == "other.cloud")] <- "Other-rare"

glmm_final_summaries_only.sig <- glmm_final_summaries[which(glmm_final_summaries$Pr...z.. < 0.05), ]

COG_category_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type == "COG category"), "variable"], decreasing = TRUE))
redundant_interaction_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type == non_redundant_category_string), "variable"], decreasing = TRUE))

glmm_final_summaries_only.sig$variable <- factor(glmm_final_summaries_only.sig$variable,
                                                 levels = rev(c("Intercept", COG_category_variables, redundant_interaction_variables)))


glmm_final_summaries$partition_clean <- factor(glmm_final_summaries$partition_clean, levels = c("Other-rare", "Shell"))


coef_barplot <- ggplot(data = glmm_final_summaries_only.sig, aes(x = Estimate, y = variable, fill = Type)) +
                        geom_bar(stat="identity") +
                        theme_bw() +
                        facet_wrap(. ~ partition_clean) +
                        ylab("Significant coefficient") +
                        geom_vline(xintercept = 0, linetype="dotted", 
                                   color = "black") +
                        theme(legend.position="bottom")

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_indepth_glmm_other.cloud_and_shell_coef.pdf',
       plot = coef_barplot,
       device = 'pdf',
       dpi = 400,
       width = 10,
       height = 10)
