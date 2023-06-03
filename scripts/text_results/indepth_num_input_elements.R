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

for (partition in c("shell", "other.cloud", "ultra.cloud")) {
  print(partition)
  element_info_subset <- element_info[which(element_info$partition == partition), ]
  print(nrow(element_info_subset))
}
