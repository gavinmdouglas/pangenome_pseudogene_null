rm(list = ls(all.names = TRUE))

# Restrict genome ids to RefSeq to get rid of some low quality assemblies
# Exception for Wolbachia pipientis too as almost all genomes for this species are Genbank-only
# And major outliers for this species are actually RefSeq ones.

tallies <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/summary/intact_and_pseudogene_tallies.tsv.gz",
                      header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE)

refseq_and_Wolbachia_pipientis <- unique(c(grep("^GCF_", rownames(tallies)),
                                           which(tallies$species == "Wolbachia_pipientis")))

tallies_mainly_refseq <- tallies[refseq_and_Wolbachia_pipientis, ]

out_ids <- data.frame(species = tallies_mainly_refseq$species, ids = rownames(tallies_mainly_refseq))

write.table(x = out_ids,
            file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/mapfiles/species_accessions_RefSeq_excluding_Wolbachia.tsv",
            sep = "\t", col.names = FALSE, row.names = FALSE, quote = FALSE)

