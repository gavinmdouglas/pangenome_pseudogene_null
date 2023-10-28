rm(list = ls(all.names = TRUE))

options(scipen = 10)

overall_p <- function(model_obj) {
  f <- summary(model_obj)$fstatistic
  p <- pf(f[1], f[2], f[3], lower.tail=FALSE)
  attributes(p) <- NULL
  return(p)
}

out_models <- readRDS(file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/model_output/pangenome_linear_models.rds')

response_variables <- c('mean_num_genes', 'genomic_fluidity', 'mean_percent_singletons_per9', 'si_sp')


model_summaries <- data.frame(matrix(NA,
                                     nrow = length(out_models$si$coefficients),
                                     ncol = 10))
colnames(model_summaries) <- c('Variable', 'Count',
                               paste(response_variables, 'estimate', sep = '_'),
                               paste(response_variables, 'P', sep = '_'))
rownames(model_summaries) <- names(out_models$si$coefficients)
model_summaries$Variable <- names(out_models$si$coefficients)

model_summaries$Variable <- gsub('class_cleanc__', '', model_summaries$Variable)
model_summaries$Variable <- gsub('dnds', 'dN/dS', model_summaries$Variable)
model_summaries$Variable <- gsub('ds', 'dS', model_summaries$Variable)
model_summaries$Variable <- gsub('_zscore', '', model_summaries$Variable)
model_summaries$Variable <- gsub('\\(Intercept\\)', 'Intercept', model_summaries$Variable)

for (response in response_variables) {
  response_summary <- summary(out_models[[response]])
  model_summaries[, paste(response, 'estimate', sep = '_')] <- format(round(response_summary$coefficients[rownames(model_summaries), 'Estimate'], digits=3), nsmall = 3)
  model_summaries[, paste(response, 'P', sep = '_')] <- response_summary$coefficients[rownames(model_summaries), 'Pr(>|t|)']
}

# Re-order dataframe rows (and indicate number of species in each class):
Intercept_row <- model_summaries['(Intercept)', ]

ds_rows <- model_summaries[grep('dS', model_summaries$Variable), ]

Class_rows <- model_summaries[grep('class', rownames(model_summaries)), ]
num_class_instances_table <- table(out_models$si$model$class_clean)
num_class_instances <- as.integer(num_class_instances_table)
names(num_class_instances) <- names(num_class_instances_table)

Class_rows[, 'Count'] <- as.integer(num_class_instances_table[gsub('class_clean', '', rownames(Class_rows))])

Class_rows <- Class_rows[order(Class_rows$Count, decreasing = TRUE), ]

Intercept_row[, 'Count'] <- as.integer(num_class_instances_table['Other'])
  
model_summaries_reordered <- do.call(rbind, list(Intercept_row, Class_rows, ds_rows))

# Add in R-squared values too.
model_summaries_reordered['rsquared', 'Variable'] <- 'Adjusted R-squared'
model_summaries_reordered['rsquared', 'Count'] <- NA
model_summaries_reordered['rsquared', paste(response_variables, 'estimate', sep = '_')] <- sapply(response_variables, function(x) { format(round(summary(out_models[[x]])$adj.r.squared, digits=3), nsmall = 3) })
model_summaries_reordered['rsquared', paste(response_variables, 'P', sep = '_')] <- sapply(response_variables, function(x) { overall_p(out_models[[x]]) })

write.table(x = model_summaries_reordered,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/model_output/linear_model_coef.tsv',
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

options(scipen = 0)
