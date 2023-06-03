rm(list = ls(all.names = TRUE))

# Spearman correlations of molecular evolution metrics vs % intact singletons

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

dnds_vs_si_cor <- cor.test(pangenome$mean_percent_singletons_per9, pangenome$median_MG94_dnds, method = "spearman")
ds_vs_si_cor <- cor.test(pangenome$mean_percent_singletons_per9, pangenome$median_pairwise_ds, method = "spearman")

round(dnds_vs_si_cor$estimate, 2)
round(dnds_vs_si_cor$p.value, 2)

round(ds_vs_si_cor$estimate, 2)
round(ds_vs_si_cor$p.value, 2)



dnds_vs_sisp_cor <- cor.test(pangenome$si_sp, pangenome$median_MG94_dnds, method = "spearman")
ds_vs_sisp_cor <- cor.test(pangenome$si_sp, pangenome$median_pairwise_ds, method = "spearman")

round(dnds_vs_sisp_cor$estimate, 2)
round(dnds_vs_sisp_cor$p.value, 2)

round(ds_vs_sisp_cor$estimate, 2)
round(ds_vs_sisp_cor$p.value, 2)