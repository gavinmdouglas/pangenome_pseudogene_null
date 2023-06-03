rm(list = ls(all.names = TRUE))

# Get % elements in significantly enriched / depleted categories. Serve as *very rough* summary statistics

library(ggplot2)

glmm_fit.info_RAW <- list()
glmm_final_summaries_RAW <- list()

# First, compare the model fits based on varying parameter complexity.
for (partition in c("shell", "other.cloud", "ultra.cloud")) {

  glmm_summaries <- readRDS(paste("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/glmm_output/", partition, "_model_output_summaries.rds", sep = ""))

  glmm_final_summaries_RAW[[partition]] <- data.frame(glmm_summaries$COG.only_species.interaction_redundant.both.interaction$coefficients$cond)
  glmm_final_summaries_RAW[[partition]]$variable <- rownames(glmm_final_summaries_RAW[[partition]])
  glmm_final_summaries_RAW[[partition]]$partition <- partition

}

glmm_final_summaries <- do.call(rbind, glmm_final_summaries_RAW)
rownames(glmm_final_summaries) <- NULL

glmm_final_summaries$Type <- NA
glmm_final_summaries[which(glmm_final_summaries$variable == "(Intercept)"), "Type"] <- "Intercept (COG category K,\nredundant COG id)"
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
glmm_final_summaries$Type[which(glmm_final_summaries$Type == "Non-redundant")] <- "Non-redundant\n+COG category:Non-redundant"

glmm_final_summaries$Type <- factor(glmm_final_summaries$Type, levels = c("Intercept (COG category K,\nredundant COG id)",
                                                                          "COG category", "Non-redundant\n+COG category:Non-redundant"))

glmm_final_summaries$partition_clean <- glmm_final_summaries$partition
glmm_final_summaries$partition_clean[which(glmm_final_summaries$partition == "shell")] <- "Shell"
glmm_final_summaries$partition_clean[which(glmm_final_summaries$partition == "other.cloud")] <- "Other-cloud"
glmm_final_summaries$partition_clean[which(glmm_final_summaries$partition == "ultra.cloud")] <- "Ultra-cloud"


glmm_final_summaries_only.sig <- glmm_final_summaries[which(glmm_final_summaries$Pr...z.. < 0.05), ]

COG_category_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type == "COG category"), "variable"], decreasing = TRUE))
redundant_interaction_variables <- sort(unique(glmm_final_summaries_only.sig[which(glmm_final_summaries_only.sig$Type =="Non-redundant\n+COG category:Non-redundant"), "variable"], decreasing = TRUE))

glmm_final_summaries_only.sig$variable <- factor(glmm_final_summaries_only.sig$variable,
                                                 levels = rev(c("Intercept", COG_category_variables, redundant_interaction_variables)))



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

# Keep only *intact* genes for this analysis.
element_info <- element_info[-which(element_info$pseudogenized), ]

COG_sig_dir <- list()
COG_sig_dir[["ultra.cloud"]] <- list()

COG_sig_dir[["ultra.cloud"]]$pos.enriched_COG_false.non.redundant <- c("C", "F", "J", "S", "X")
COG_sig_dir[["ultra.cloud"]]$neg.enriched_COG_false.non.redundant <- c("D")

COG_sig_dir[["ultra.cloud"]]$pos.enriched_COG_true.non.redundant <- as.character()
COG_sig_dir[["ultra.cloud"]]$neg.enriched_COG_true.non.redundant <- c("E", "F", "H", "J", "L", "N", "P", "Q" ,"R", "S", "U", "V", "W", "X")


COG_sig_dir[["other.cloud"]] <- list()

COG_sig_dir[["other.cloud"]]$pos.enriched_COG_false.non.redundant <- c("X")
COG_sig_dir[["other.cloud"]]$neg.enriched_COG_false.non.redundant <- as.character()

COG_sig_dir[["other.cloud"]]$pos.enriched_COG_true.non.redundant <- "C"
COG_sig_dir[["other.cloud"]]$neg.enriched_COG_true.non.redundant <- c("D", "E", "F", "G", "H", "I", "J", "L", "M", "N", "O", "P", "Q" ,"R", "S", "T", "U", "V", "W", "X")


COG_sig_dir[["shell"]] <- list()

COG_sig_dir[["shell"]]$pos.enriched_COG_false.non.redundant <- as.character()
COG_sig_dir[["shell"]]$neg.enriched_COG_false.non.redundant <- c("I", "N", "W")

COG_sig_dir[["shell"]]$pos.enriched_COG_true.non.redundant <- c("S", "V", "W")
COG_sig_dir[["shell"]]$neg.enriched_COG_true.non.redundant <- c("C", "E", "F", "G", "H", "I", "J", "L", "M", "N", "O", "P", "Q", "R", "T", "U", "X")

species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene/mapfiles/species.txt", sep = "\t")$V1
summary_tab_RAW <- list()
for (partition in names(COG_sig_dir)) {
  
  summary_tab_RAW[[partition]] <- data.frame(matrix(NA, nrow = 10, ncol = 3))
  rownames(summary_tab_RAW[[partition]]) <- species
  colnames(summary_tab_RAW[[partition]]) <- c("pseudo_enriched", "pseudo_depleted", "other")
  
  element_info_subset <- element_info[which(element_info$partition == partition), ]
   
  for (sp in species) {
    
    element_info_subset_sp <- element_info_subset[which(element_info_subset$species == sp), , drop = FALSE]
    
    total_elements <- nrow(element_info_subset_sp)
    
    summary_tab_RAW[[partition]][sp, "pseudo_enriched"] <- nrow(element_info_subset_sp[which(element_info_subset_sp$COG_category %in% COG_sig_dir[[partition]]$pos.enriched_COG_false.non.redundant & ! element_info_subset_sp$NO_redundant_intact_COG), ]) +
                                                           nrow(element_info_subset_sp[which(element_info_subset_sp$COG_category %in% COG_sig_dir[[partition]]$pos.enriched_COG_true.non.redundant & element_info_subset_sp$NO_redundant_intact_COG), ])
    
    summary_tab_RAW[[partition]][sp, "pseudo_depleted"] <- nrow(element_info_subset_sp[which(element_info_subset_sp$COG_category %in% COG_sig_dir[[partition]]$neg.enriched_COG_false.non.redundant & ! element_info_subset_sp$NO_redundant_intact_COG), ]) +
                                                           nrow(element_info_subset_sp[which(element_info_subset_sp$COG_category %in% COG_sig_dir[[partition]]$neg.enriched_COG_true.non.redundant & element_info_subset_sp$NO_redundant_intact_COG), ])

    summary_tab_RAW[[partition]][sp, "other"] <- total_elements - summary_tab_RAW[[partition]][sp, "pseudo_enriched"] - summary_tab_RAW[[partition]][sp, "pseudo_depleted"]
  }
  
}

summary_tab_percent <- list()

for (partition in names(summary_tab_RAW)) {
  summary_tab_percent[[partition]] <- data.frame(sweep(x = summary_tab_RAW[[partition]], MARGIN = 1, STATS = rowSums(summary_tab_RAW[[partition]]), FUN = '/')) * 100
}

round(colMeans(summary_tab_percent$ultra.cloud), 2)
round(colMeans(summary_tab_percent$other.cloud), 2)
round(colMeans(summary_tab_percent$shell), 2)

round(sapply(summary_tab_percent$ultra.cloud, sd), 2)
round(sapply(summary_tab_percent$other.cloud, sd), 2)
round(sapply(summary_tab_percent$shell, sd), 2)


# Also, to get the % of non-redundant elements across all species:
round((sum(element_info$NO_redundant_intact_COG) / nrow(element_info)) * 100, 2)

