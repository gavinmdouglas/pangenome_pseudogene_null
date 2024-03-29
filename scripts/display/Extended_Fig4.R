rm(list = ls(all.names = TRUE))

library(ggplot2)

COG_to_category <- read.table('/data1/gdouglas/db/COG_definitions/COG_category_descrip_succinct.tsv',
                              header = FALSE, sep = '\t', stringsAsFactors = FALSE, row.names = 1)

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
  
  glmm_final_summaries[which(glmm_final_summaries$partition == partition_level & glmm_final_summaries$Type == "Non-redundant"), "Std..Error"] <-
    glmm_final_summaries[which(glmm_final_summaries$partition == partition_level & glmm_final_summaries$Type == "Non-redundant"), "Std..Error"] +
    glmm_final_summaries[which(glmm_final_summaries$partition == partition_level & glmm_final_summaries$variable == "Non-redundant"), "Std..Error"]
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
orig_COG_category_variables <- COG_category_variables
COG_category_variables <- paste(COG_category_variables, '-', COG_to_category[COG_category_variables, 'V2'])

redundant_interaction_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type == non_redundant_category_string), "variable"], decreasing = TRUE))
orig_redundant_interaction_variables <- redundant_interaction_variables

redundant_interaction_variables <- gsub(' \\(Non-redundant\\)', '', redundant_interaction_variables)
redundant_interaction_variables <- paste(redundant_interaction_variables, '-', COG_to_category[redundant_interaction_variables, 'V2'], '(Non-redun.)')

variable_rename <- data.frame(old = c('Intercept', orig_COG_category_variables, orig_redundant_interaction_variables),
                              new = c('Intercept', COG_category_variables, redundant_interaction_variables))
rownames(variable_rename) <- variable_rename$old

for (i in 1:nrow(glmm_final_summaries_only.sig)) {
  glmm_final_summaries_only.sig[i, 'variable'] <- variable_rename[glmm_final_summaries_only.sig[i, 'variable'], 'new']
}


glmm_final_summaries_only.sig$variable <- factor(glmm_final_summaries_only.sig$variable,
                                                 levels = rev(c("Intercept", COG_category_variables, redundant_interaction_variables)))


glmm_final_summaries$partition_clean <- factor(glmm_final_summaries$partition_clean, levels = c("Other-rare", "Shell"))


coef_barplot <- ggplot(data = glmm_final_summaries_only.sig, aes(x = Estimate, y = variable, fill = Type)) +
                        geom_bar(stat="identity") +
                        theme_bw() +
                        scale_fill_manual(values = c("#b16450", "#a0c68f", "#4f5f42")) +
                        facet_wrap(. ~ partition_clean) +
                        ylab("Significant coefficient") +
                        geom_vline(xintercept = 0, linetype="dotted", 
                                   color = "black") +
                        theme(legend.position="bottom") +
                        geom_errorbar(aes(xmin = Estimate - Std..Error,
                                          xmax = Estimate + Std..Error),
                                      width = 0.2, color = "black")

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/Douglas_ED_Fig4.pdf',
       plot = coef_barplot,
       device = 'pdf',
       dpi = 400,
       width = 10,
       height = 10)


# Write out source data:
source_out <- glmm_final_summaries_only.sig[, c('partition', 'variable', 'Type', 'Estimate', 'Std..Error', 'Pr...z..')]
colnames(source_out) <- c('partition', 'variable', 'type', 'coefficient', 'standard_error', 'P_value')
source_out$type <- as.character(source_out$type)
source_out$type <- gsub('\n', ' ', source_out$type)
write.table(x = source_out,
            file = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data/ED_Fig4.tsv",
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
