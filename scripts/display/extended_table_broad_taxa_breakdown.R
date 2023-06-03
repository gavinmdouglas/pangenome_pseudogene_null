rm(list = ls(all.names = TRUE))

# Info on genomes for broad analysis.
pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

sp_taxonomy <- read.table('/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/taxonomy.tsv.gz',
                          row.names = 1, stringsAsFactors = FALSE, sep = '\t', header = TRUE)
sp_taxonomy <- sp_taxonomy[rownames(pangenome), ]

# Generate a table with quick breakdowns of taxa counts at various high levels.
prep_character_tally_line <- function(tallies) {
  tallies <- sort(tallies, decreasing = TRUE)
  working <- as.character()
  for (i in 1:length(tallies)) {
    tally_name <- names(tallies)[i]
    tally <- tallies[i]
    working <- c(working, paste(tally_name, ': ', as.character(tally), sep = ''))
  }
  return(paste(working, collapse = '; ', sep = ''))
}

domain_counts <- c(length(grep("Bacteria", sp_taxonomy$gtdb_taxonomy)),
                   length(grep("Archaea", sp_taxonomy$gtdb_taxonomy)))
names(domain_counts) <- c("Bacteria", "Archaea")

phylum_counts <- table(sp_taxonomy$phylum)
class_counts <- table(sp_taxonomy$class)

names(phylum_counts) <- gsub('^p__', '', names(phylum_counts))
names(class_counts) <- gsub('^c__', '', names(class_counts))

taxa_breakdown <- data.frame(Level = c('Domain', "Phylum", 'Class'),
                             Breakdown = c(prep_character_tally_line(domain_counts),
                                           prep_character_tally_line(phylum_counts),
                                           prep_character_tally_line(class_counts)))

write.table(x = taxa_breakdown,
            file = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_table_broad_taxa_breakdown.tsv',
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
