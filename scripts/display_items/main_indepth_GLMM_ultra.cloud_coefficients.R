rm(list = ls(all.names = TRUE))

library(ggplot2)

glmm_fit.info_RAW <- list()
glmm_final_summaries_RAW <- list()

# Parse coefficents for ultra-cloud only.
partition <- "ultra.cloud"

glmm_summaries <- readRDS(paste("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/glmm_output/", partition, "_model_output_summaries.rds", sep = ""))

glmm_final_summaries_RAW[[partition]] <- data.frame(glmm_summaries$COG.only_species.interaction_redundant.both.interaction$coefficients$cond)
glmm_final_summaries_RAW[[partition]]$variable <- rownames(glmm_final_summaries_RAW[[partition]])
glmm_final_summaries_RAW[[partition]]$partition <- partition

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

glmm_final_summaries_only.sig <- glmm_final_summaries[which(glmm_final_summaries$Pr...z.. < 0.05), ]

COG_category_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type == "COG category"), "variable"], decreasing = TRUE))
redundant_interaction_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type == non_redundant_category_string), "variable"], decreasing = TRUE))

glmm_final_summaries_only.sig$variable <- factor(glmm_final_summaries_only.sig$variable,
                                                 levels = rev(c("Intercept", COG_category_variables, redundant_interaction_variables)))


ultra.cloud_coef_barplot <- ggplot(data = glmm_final_summaries_only.sig, aes(x = Estimate, y = variable, fill = Type)) +
                                    geom_bar(stat="identity") +
                                    theme_bw() +
                                    ylab("Significant coefficient") +
                                    geom_vline(xintercept = 0, linetype="dotted", 
                                               color = "black")

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_indepth_glmm_ultra.cloud_coef.pdf',
       plot = ultra.cloud_coef_barplot,
       device = 'pdf',
       dpi = 600,
       width = 8,
       height = 6)
