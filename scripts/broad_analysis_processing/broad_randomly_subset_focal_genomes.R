rm(list = ls(all.names = TRUE))

cluster_member_breakdown <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene_PRIOR/summary/cluster_member_breakdown.tsv.gz",
                                       header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 2)

species <- read.table("/data1/gdouglas/projects/accessory_vs_pseudogene_PRIOR/mapfiles/species.txt",
                      stringsAsFactors = FALSE)$V1

for (sp in species) {
 
  sp_accessions <- unique(cluster_member_breakdown[which(cluster_member_breakdown$species == sp), "accession"])
   
  sp_random_accessions <- sample(unique(sp_accessions), 20)
  
  outfile <- paste('/data1/gdouglas/projects/accessory_vs_pseudogene/Ne_based_analysis/focal_genomes_subsampled/subsampled_genome_ids/',
                   sp, '.txt', sep = '')
  
  write.table(x = sp_random_accessions, file = outfile, row.names = FALSE, col.names = FALSE, quote = FALSE)
  
}