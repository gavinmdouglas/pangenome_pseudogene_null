rm(list = ls(all.names = TRUE))

# For non-focal genomes.

outfiles <- list.files("/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/panaroo_non.focal_cds/pairwise_dnds_stats/", full.names = TRUE)

species <-  gsub(".tsv.gz", "", list.files("/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/panaroo_non.focal_cds/pairwise_dnds_stats/"))

dn <- as.numeric()
ds <- as.numeric()
dnds <- as.numeric()

for (i in 1:length(outfiles)) {
  sp_out <- read.table(file = outfiles[i], header = FALSE, stringsAsFactors = FALSE, row.names = 1, sep = "\t")
  colnames(sp_out) <- c('mean_n_subs', 'mean_n_sites', 'mean_s_subs', 'mean_s_sites', 'mean_dn', 'mean_ds', 'mean_dnds')

  dn <- c(dn, median(sp_out$mean_dn))
  ds <- c(ds, median(sp_out$mean_ds))
  dnds <- c(dnds, median(sp_out$mean_dnds, na.rm = TRUE))  
  
}

output <- data.frame(species = species, dn = dn, ds = ds, dnds = dnds)

write.table(x = output,
            file = "/data1/gdouglas/projects/accessory_vs_pseudogene/non.focal_genomes/panaroo_non.focal_cds/median_pairwise_dnds_stats.tsv",
            sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
