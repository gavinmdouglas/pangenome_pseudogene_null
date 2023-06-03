rm(list = ls(all.names = TRUE))

# Heatmaps summarizing percentages + counts of each species' pangenome found in each partition of the pangenome.
# Also show breakdown of % clusters without COG category annotation per cell.

library(ComplexHeatmap)
library(circlize)
library(cowplot)
library(ggplot2)

# % of clusters by pangenome partition (based on *all* clusters)
all.clusters_percent_by_type <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_breakdown_tables/all.clusters_percent_by_type.tsv.gz",
                                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

all.clusters_percent_by_type_fill <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_breakdown_tables/all.clusters_percent_by_type_fill.tsv.gz",
                                                      header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

all.clusters_percent_by_type_fill <- apply(X = all.clusters_percent_by_type_fill,
                                                 MARGIN = 2,
                                                 FUN = function(x) { gsub(" \\(", "\n\\(", x) })

all.clusters_percent_by_type_heatmap <- Heatmap(matrix = as.matrix(all.clusters_percent_by_type),
                                                      name = "Percent",
                                                      row_names_side = "left",
                                                      row_labels = gt_render(paste("*", gsub("_", " ", rownames(all.clusters_percent_by_type)), "*", sep = "")),
                                                      column_labels = c("Ultra-rare", "Other-rare", "Shell", "Soft-core", "Ultra-rare", "Other-rare", "Shell", "Soft-core", "Ultra-rare", "Other-rare", "Shell", "Soft-core"),
                                                      column_split = c("Intact", "Intact", "Intact", "Intact", "Mixed", "Mixed", "Mixed", "Mixed", "Pseudogene", "Pseudogene", "Pseudogene", "Pseudogene"),
                                                      cluster_rows = FALSE,
                                                      column_gap = unit(5, "mm"),
                                                      cluster_columns = FALSE,
                                                      column_names_rot = 45,
                                                      col = colorRamp2(c(0, 50, 100), c("slateblue1","white", "red")),
                                                      cell_fun = function(j, i, x, y, width, height, fill) {
                                                        grid.text(all.clusters_percent_by_type_fill[i, j], x, y, gp = gpar(fontsize = 10, fontface="bold"))
                                                      })

# % of clusters in above heatmap that are *not* COG-annotated.
unannot_percent_by_type.cell <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_breakdown_tables/unannot_percent_by_table.cell.tsv.gz",
                                           header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

unannot_percent_by_type.cell_fill <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_breakdown_tables/unannot_percent_by_table.cell_fill.tsv.gz",
                                                header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

unannot_percent_by_type.cell_fill <- apply(X = unannot_percent_by_type.cell_fill,
                                           MARGIN = 2,
                                           FUN = function(x) { gsub(" \\(", "\n\\(", x) })

unannot_percent_by_type.cell_heatmap <- Heatmap(matrix = as.matrix(unannot_percent_by_type.cell),
                                                name = "Percent",
                                                row_names_side = "left",
                                                row_labels = gt_render(paste("*", gsub("_", " ", rownames(unannot_percent_by_type.cell)), "*", sep = "")),
                                                column_labels = c("Ultra-rare", "Other-rare", "Shell", "Soft-core", "Ultra-rare", "Other-rare", "Shell", "Soft-core", "Ultra-rare", "Other-rare", "Shell", "Soft-core"),
                                                column_split = c("Intact", "Intact", "Intact", "Intact", "Mixed", "Mixed", "Mixed", "Mixed", "Pseudogene", "Pseudogene", "Pseudogene", "Pseudogene"),
                                                cluster_rows = FALSE,
                                                column_gap = unit(5, "mm"),
                                                cluster_columns = FALSE,
                                                column_names_rot = 45,
                                                col = colorRamp2(c(0, 50, 100), c("slateblue1","white", "red")),
                                                cell_fun = function(j, i, x, y, width, height, fill) {
                                                  grid.text(unannot_percent_by_type.cell_fill[i, j], x, y, gp = gpar(fontsize = 10, fontface="bold"))
                                                })

# % of clusters by pangenome partition (based on COG-category-annotated clusters only)
cog.category.annot_percent_by_type <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_breakdown_tables/cog.category.annot_percent_by_type.tsv.gz",
                                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

cog.category.annot_percent_by_type_fill <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/cluster_breakdown_tables/cog.category.annot_percent_by_type_fill.tsv.gz",
                                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

cog.category.annot_percent_by_type_fill <- apply(X = cog.category.annot_percent_by_type_fill,
                                                 MARGIN = 2,
                                                 FUN = function(x) { gsub(" \\(", "\n\\(", x) })

cog.category.annot_percent_by_type_heatmap <- Heatmap(matrix = as.matrix(cog.category.annot_percent_by_type),
                                                      name = "Percent",
                                                      row_names_side = "left",
                                                      row_labels = gt_render(paste("*", gsub("_", " ", rownames(cog.category.annot_percent_by_type)), "*", sep = "")),
                                                      column_labels = c("Ultra-rare", "Other-rare", "Shell", "Soft-core", "Ultra-rare", "Other-rare", "Shell", "Soft-core", "Ultra-rare", "Other-rare", "Shell", "Soft-core"),
                                                      column_split = c("Intact", "Intact", "Intact", "Intact", "Mixed", "Mixed", "Mixed", "Mixed", "Pseudogene", "Pseudogene", "Pseudogene", "Pseudogene"),
                                                      cluster_rows = FALSE,
                                                      column_gap = unit(5, "mm"),
                                                      cluster_columns = FALSE,
                                                      column_names_rot = 45,
                                                      col = colorRamp2(c(0, 50, 100), c("slateblue1","white", "red")),
                                                      cell_fun = function(j, i, x, y, width, height, fill) {
                                                        grid.text(cog.category.annot_percent_by_type_fill[i, j], x, y, gp = gpar(fontsize = 10, fontface="bold"))
                                                      })

all.cluster_and_unannot_heatmaps <- plot_grid(grid.grabExpr(draw(all.clusters_percent_by_type_heatmap,
                                                                 column_title = 'Percent of all clusters (including those lacking COG category annotation) per pangenome partition')),
                                              grid.grabExpr(draw(unannot_percent_by_type.cell_heatmap,
                                                                 column_title = 'Percent of clusters in each cell of the above panel that lack a COG category annotation')),
                                              labels = c('a', 'b'),
                                              nrow = 2)

cog.category.annot_percent_by_type_heatmap <- plot_grid(grid.grabExpr(draw(cog.category.annot_percent_by_type_heatmap,
                                                                           column_title = 'Percent of clusters (restricted to only those with COG category annotations) per pangenome partition')))

ggsave(plot = all.cluster_and_unannot_heatmaps,
       filename = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_indepth_all.cluster_and_unannot_heatmaps.pdf",
       device = "pdf", width = 12, height = 12, units = "in", dpi = 400)

ggsave(plot = cog.category.annot_percent_by_type_heatmap,
       filename = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_indepth_cog.category.annot_percent_by_type_heatmap.pdf",
       device = "pdf", width = 12, height = 6, units = "in", dpi = 400)
