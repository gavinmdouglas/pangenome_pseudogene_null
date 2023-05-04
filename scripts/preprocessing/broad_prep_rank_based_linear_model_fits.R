rm(list = ls(all.names = TRUE))

options(scipen = 10)

out_models <- readRDS(file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/model_output/pangenome_rank_based_linear_models.rds')

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
  model_summaries[, paste(response, 'P', sep = '_')] <- format(round(response_summary$coefficients[rownames(model_summaries), 'p.value'], digits=3), nsmall = 3)
}

# Re-order dataframe rows (and indicate number of species in each class):
Intercept_row <- model_summaries['(Intercept)', ]
ds_rows <- model_summaries[grep('dS', model_summaries$Variable), ]

Class_rows <- model_summaries[grep('class', rownames(model_summaries)), ]
num_class_instances <- colSums(out_models$si$x[, grep('class', colnames(out_models$si$x))])
Class_rows[names(num_class_instances), 'Count'] <- as.integer(num_class_instances)
Class_rows <- Class_rows[order(Class_rows$Count, decreasing = TRUE), ]

model_summaries_reordered <- do.call(rbind, list(Intercept_row, Class_rows, ds_rows))

write.table(x = model_summaries_reordered,
            file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/model_output/rank_based_linear_coef.tsv',
            sep = '\t', row.names = FALSE, col.names = TRUE, quote = FALSE)

options(scipen = 0)
