rm(list = ls(all.names = TRUE))

# Fit linear models on pangenome diversity metrics vs molecular evolution metrics and taxonomic classes.
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

pangenome <- pangenome[which(pangenome$mean_percent_singletons_per9 > 0), ]

# Classes with <= 5 members were collapsed into "Other". These genomes were used as the first level of the class variable, i.e., used for the intercept.
class_tallies <- table(pangenome$class)
rare_classes <- names(class_tallies)[which(class_tallies <= 5)]
pangenome$class_clean <- pangenome$class
pangenome$class_clean[which(pangenome$class %in% rare_classes)] <- "Other"
class_order <- sort(unique(pangenome$class_clean))
pangenome$class_clean <- factor(pangenome$class_clean, levels = c("Other", class_order[-which(class_order == "Other")]))

# Mean centre and scale dS and dN/dS
pangenome$dnds <- -1 * (1 / pangenome$dnds)
pangenome$ds <- sqrt(pangenome$ds)
pangenome$dnds_zscore <- (pangenome$dnds - mean(pangenome$dnds)) / sd(pangenome$dnds)
pangenome$ds_zscore <- (pangenome$ds - mean(pangenome$ds)) / sd(pangenome$ds)

out_models <- list()

response_variables <- c('mean_num_genes', 'genomic_fluidity', 'mean_percent_singletons_per9', 'si_sp')

for (response in response_variables) {

  # Mean no. genes already roughly normal, and becomes bimodal when transformed, so leave it unchanged. 
  if (response != 'mean_num_genes') {
    pangenome[, response] <- sqrt(pangenome[, response])
  }

  pangenome[, response] <- (pangenome[, response] - mean(pangenome[, response])) / sd(pangenome[, response])

  out_models[[response]] <- lm(formula = paste0(response, ' ~ class_clean + ds_zscore + dnds_zscore'), data = pangenome)

}

saveRDS(object = out_models,
        file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/model_output/pangenome_linear_models.rds')
