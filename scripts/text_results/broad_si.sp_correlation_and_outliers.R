rm(list = ls(all.names = TRUE))

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Pearson correlation
si_vs_sp_pearson <- cor.test(pangenome$mean_percent_singletons_per9,
                             pangenome$mean_percent_singletons_pseudo_per9,
                             method = "pearson")

round(si_vs_sp_pearson$estimate, 2)
si_vs_sp_pearson$p.value

# Species with max and min ratios, mentioned in text:
pangenome[which.max(pangenome$si_sp), ]
pangenome[which.min(pangenome$si_sp), ]
