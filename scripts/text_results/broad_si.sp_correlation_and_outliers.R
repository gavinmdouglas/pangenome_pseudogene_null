rm(list = ls(all.names = TRUE))

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Spearman correlation
si_vs_sp_spearman <- cor.test(pangenome$mean_percent_singletons_per9,
                             pangenome$mean_percent_singletons_pseudo_per9,
                             method = "spearman")

round(si_vs_sp_spearman$estimate, 2)
si_vs_sp_spearman$p.value

# Species with max and min ratios, mentioned in text:
pangenome[which.max(pangenome$si_sp), ]
pangenome[which.min(pangenome$si_sp), ]
