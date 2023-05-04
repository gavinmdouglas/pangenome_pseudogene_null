rm(list = ls(all.names = TRUE))

library(ComplexHeatmap)

model_summaries <- read.table(file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/model_output/linear_model_coef.tsv.gz',
                              sep = '\t', stringsAsFactors = FALSE, header = TRUE)

estimates <- model_summaries[, grep('estimate', colnames(model_summaries))]

pvalues <- model_summaries[, grep('_P$', colnames(model_summaries))]
pvalues_discrete <- pvalues
pvalues_discrete[pvalues < 0.05] <-  '< 0.05'
pvalues_discrete[pvalues >= 0.05] <- '>= 0.05'

clean_variables <- colnames(estimates)
clean_variables[which(clean_variables == 'mean_num_genes_estimate')] <- 'Mean no. genes'
clean_variables[which(clean_variables == 'genomic_fluidity_estimate')] <- 'Genomic fluidity'
clean_variables[which(clean_variables == 'mean_percent_singletons_per9_estimate')] <- expression(s[i])
clean_variables[which(clean_variables == 'si_sp_estimate')] <- expression(s[i]/s[p])

class_count <- model_summaries$Count

class_count_annotation <- rowAnnotation(`No.\ngenomes\nin class` = class_count,
                                        col = list(`No.\ngenomes\nin class` = circlize::colorRamp2(c(0,
                                                                                                     min(class_count, na.rm = TRUE),
                                                                                                     max(class_count, na.rm = TRUE)),
                                                                                                   c("grey", "skyblue1", "skyblue4"))),
                                        na_col = 'white',
                                        show_annotation_name = FALSE)

model_summaries$Variable[which(model_summaries$Variable == "Adjusted R-squared")] <- "\"Adjusted R\"^2"

row_splits <- factor(c('Intercept\ncoef.', 'Intercept\ncoef.', 'Intercept\ncoef.', 'Intercept\ncoef.',
                       'Intercept\ncoef.', 'Intercept\ncoef.', 'Intercept\ncoef.', 'Intercept\ncoef.',
                       'Slope\ncoef.', 'Slope\ncoef.',
                       'Model\nfit'),
                     levels = c('Intercept\ncoef.', 'Slope\ncoef.', 'Model\nfit'))

heatmap_linear_coef_raw <- Heatmap(matrix = as.matrix(pvalues_discrete),
                                          name = 'P-value',
                                          row_names_side = "left",
                                          column_labels = clean_variables,
                                          row_labels = parse(text = model_summaries$Variable),
                                          row_split = row_splits,
                                          cluster_columns = FALSE,
                                          cluster_rows = FALSE,
                                          cluster_row_slices = FALSE,
                                          column_names_rot = 45,
                                          row_title_rot = 0,
                                          row_gap=unit(0.05, "npc"),
                                          right_annotation = class_count_annotation,
                                          col = c("darkseagreen2", "gray95"),
                                          cell_fun = function(j, i, x, y, width, height, fill) {
                                            grid.text(format(round(estimates[i, j], digits=3), nsmall = 3), x, y, gp = gpar(fontsize = 10, fontface="bold"))
                                          })

heatmap_linear_coef <- cowplot::plot_grid(grid.grabExpr(draw(heatmap_linear_coef_raw)))

ggsave(plot = heatmap_linear_coef,
       filename = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/main_linear_model_coef.pdf",
       device = "pdf", width = 6, height = 5, units = "in", dpi = 600)
