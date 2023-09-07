rm(list = ls(all.names = TRUE))

library(cowplot)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)

run_pairwise_pangenome_cor <- function(metrics, in_df) {
  
  spearman_rho <- data.frame(matrix(NA,
                                    ncol = length(metrics),
                                    nrow = length(metrics)))
  rownames(spearman_rho) <- metrics
  colnames(spearman_rho) <- metrics
  spearman_p <- spearman_rho
  
  for (i in 1:(length(metrics) - 1)) {
    
    if (i == length(metrics)) { next }
    
    for (j in (i + 1):length(metrics)) {
      
      if (metrics[j] == 'dnds_by_ds') { next }
      
      metric_j_vals <- in_df[, metrics[j]]
      
      if (metrics[i] == 'dnds_by_ds') {
        
        partial_spearman_out <- ppcor::pcor.test(x = in_df$dnds,
                                                 y = metric_j_vals,
                                                 z = in_df$ds,
                                                 method = 'spearman')
        
        spearman_rho[metrics[j], 'dnds_by_ds'] <- partial_spearman_out$estimate
        
        spearman_p[metrics[j], 'dnds_by_ds'] <- partial_spearman_out$p.value
        
      } else {
        
        metric_i_vals <- in_df[, metrics[i]]
        
        spearman_out <- cor.test(metric_i_vals,
                                 metric_j_vals,
                                 method = 'spearman',
                                 exact = FALSE)
        
        spearman_rho[metrics[j], metrics[i]] <- spearman_out$estimate
        
        spearman_p[metrics[j], metrics[i]] <- spearman_out$p.value
        
      }
      
    }
    
  }
  
  spearman_rho_sig_only <- spearman_rho
  spearman_rho_sig_only[spearman_p >= 0.05] <- NA
  
  return(list(spearman_rho_sig_only = spearman_rho_sig_only,
              spearman_rho = spearman_rho,
              spearman_p = spearman_p))
  
}

pangenome <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz',
                        header = TRUE, sep = '\t', row.names = 1)

metrics <- c('ds',
             'dnds',
             'dnds_by_ds',
             'mean_num_genes',
             'genomic_fluidity',
             'mean_percent_singletons_per9',
             'mean_percent_singletons_pseudo_per9',
             'si_sp')

full_data_cor <- run_pairwise_pangenome_cor(metrics, pangenome)

pangenome_no.outlier <- pangenome
pangenome_no.outlier <- pangenome_no.outlier[which(pangenome_no.outlier$dnds < 0.5), ]
pangenome_no.outlier <- pangenome_no.outlier[which(pangenome_no.outlier$ds > 0.01), ]

no.outlier_data_cor <- run_pairwise_pangenome_cor(metrics, pangenome_no.outlier)


# Remove blank rows and columns
full_data_cor$spearman_rho_sig_only <- full_data_cor$spearman_rho_sig_only[-c(1, 3), 1:(length(metrics) - 1)]
full_data_cor$spearman_rho <- full_data_cor$spearman_rho[-c(1, 3), 1:(length(metrics) - 1)]
full_data_cor$spearman_p <- full_data_cor$spearman_p[-c(1, 3), 1:(length(metrics) - 1)]

no.outlier_data_cor$spearman_rho_sig_only <- no.outlier_data_cor$spearman_rho_sig_only[-c(1, 3), 1:(length(metrics) - 1)]
no.outlier_data_cor$spearman_rho <- no.outlier_data_cor$spearman_rho[-c(1, 3), 1:(length(metrics) - 1)]
no.outlier_data_cor$spearman_p <- no.outlier_data_cor$spearman_p[-c(1, 3), 1:(length(metrics) - 1)]


final_label_map <- list()
final_label_map[['ds']] <- 'dS'
final_label_map[['dnds']] <- 'dN/dS'
final_label_map[['dnds_by_ds']] <- 'dN/dS (dS partial)'
final_label_map[['mean_num_genes']] <- 'Mean no. genes'
final_label_map[['genomic_fluidity']] <- 'Genomic fluidity'
final_label_map[['mean_percent_singletons_per9']] <- 's<sub>i</sub>'
final_label_map[['mean_percent_singletons_pseudo_per9']] <- 's<sub>p</sub>'
final_label_map[['si_sp']] <- 's<sub>i</sub>/s<sub>p</sub>'

clean_row_labels <- sapply(rownames(full_data_cor$spearman_rho_sig_only), function(x) { return(final_label_map[[x]]) })
clean_column_labels <- sapply(colnames(full_data_cor$spearman_rho_sig_only), function(x) { return(final_label_map[[x]]) })

all_data_heatmap <- Heatmap(matrix = as.matrix(full_data_cor$spearman_rho_sig_only),
                            #name = 'Spearman\'s\ncorrelation\ncoefficient\n(sig. only)', 
                            show_heatmap_legend = FALSE,
                            col = colorRamp2(c(-1, 0, 1), c("tomato1","white", "slateblue1")),
                            na_col = 'gray50',
                            rect_gp = gpar(type = "none"),
                            cluster_rows = FALSE,
                            cluster_columns = FALSE,
                            column_names_rot = 45,
                            row_title_rot = 0,
                            column_title='Full dataset (668 species)',
                            row_labels = ComplexHeatmap::gt_render(clean_row_labels),
                            column_labels = ComplexHeatmap::gt_render(clean_column_labels),
                            cell_fun = function(j, i, x, y, w, h, fill) {
                              if(! is.na(full_data_cor$spearman_rho[i, j])) {
                                grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                                grid.text(sprintf("%.2f", full_data_cor$spearman_rho[i, j]), x, y, gp = gpar(fontsize = 10))
                              }
                            })


no.outlier_data_heatmap <- Heatmap(matrix = as.matrix(no.outlier_data_cor$spearman_rho_sig_only),
                                    name = 'Significant\nSpearman\'s\ncorrelation\ncoefficient',
                                    col = colorRamp2(c(-1, 0, 1), c("tomato1","white", "slateblue1")),
                                    na_col = 'gray50',
                                    rect_gp = gpar(type = "none"),
                                    cluster_rows = FALSE,
                                    cluster_columns = FALSE,
                                    column_names_rot = 45,
                                    row_title_rot = 0,
                                    column_title='Outlier species (dS <= 0.01 and\ndN/dS >= 0.5) removed (596 species)',
                                    row_labels = ComplexHeatmap::gt_render(clean_row_labels),
                                    column_labels = ComplexHeatmap::gt_render(clean_column_labels),
                                    cell_fun = function(j, i, x, y, w, h, fill) {
                                      if(! is.na(no.outlier_data_cor$spearman_rho[i, j])) {
                                        grid.rect(x, y, w, h, gp = gpar(fill = fill, col = fill))
                                        grid.text(sprintf("%.2f", no.outlier_data_cor$spearman_rho[i, j]), x, y, gp = gpar(fontsize = 10))
                                      }
                                    })

combined <- plot_grid(grid.grabExpr(draw(all_data_heatmap + no.outlier_data_heatmap,
                                         ht_gap = unit(20, "mm"))))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/WORKING_all_pairwise_cor.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 600,
       width = 12,
       height = 6)