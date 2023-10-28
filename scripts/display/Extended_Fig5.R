# Summarize OR's of individual COGs within X (mobilome) COG category based on what type of MGE is associated.

rm(list = ls(all.names = TRUE))

library(ggplot2)
library(ggbeeswarm)

COG_id_descrip <- read.table("/data1/gdouglas/db/COG_definitions/cog-20.def.tab",
                             header = FALSE, sep = "\t", stringsAsFactors = FALSE, quote = "", comment.char = "", row.names = 1)

mobilome_COG_hits <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/indepth_10_species_analysis/COG_enrichment_results/ultra.cloud-COG-gene-enrichments.tsv.gz",
                                   header = TRUE, sep = "\t", stringsAsFactors = FALSE)
mobilome_COG_hits <- mobilome_COG_hits[which(mobilome_COG_hits$fdr < 0.05), ]
mobilome_COG_hits <- mobilome_COG_hits[which(mobilome_COG_hits$COG_category == "X"), ]

# Not that this is restricted to the redundant set.
mobilome_COG_hits <- mobilome_COG_hits[which(mobilome_COG_hits$type == "pos.enriched_COG_false.non.redundant"), ]

mobilome_COG_hits$Description <- tolower(COG_id_descrip[mobilome_COG_hits$category, "V3"])

transposon_hits <- grep("transpos", mobilome_COG_hits$Description)
phage_hits <- grep("phage", mobilome_COG_hits$Description)
plasmid_hits <- unique(c(grep("plasmid", mobilome_COG_hits$Description), grep("conju", mobilome_COG_hits$Description)))

mobilome_COG_hits$Subtype <- "Mixed/other"
mobilome_COG_hits[transposon_hits[which(! transposon_hits %in% c(phage_hits, plasmid_hits))], "Subtype"] <- "Transposon"
mobilome_COG_hits[phage_hits[which(! phage_hits %in% c(transposon_hits, plasmid_hits))], "Subtype"] <- "Phage"
mobilome_COG_hits[plasmid_hits[which(! plasmid_hits %in% c(transposon_hits, phage_hits))], "Subtype"] <- "Plasmid"

mobilome_COG_hits$log2OR <- log2(mobilome_COG_hits$OR + 0.1)

# Transposon summary statistics to report in text:
mean(mobilome_COG_hits[which(mobilome_COG_hits$Subtype == 'Transposon'), 'OR'])
sd(mobilome_COG_hits[which(mobilome_COG_hits$Subtype == 'Transposon'), 'OR'])

# Sample sizes for figure legend:
table(mobilome_COG_hits$Subtype)

mobilome_enrichment_by_subtype <- ggplot(data = mobilome_COG_hits, aes(x = Subtype, y = log2OR)) +
                                        geom_boxplot(outlier.shape = NA) +
                                        geom_quasirandom(alpha = 0.8) +
                                        theme_bw() +
                                        xlab("Associated mobile genetic element sub-type") +
                                        ylab(expression('log'[2]*'(Odd\'s ratio + 0.1)'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/Douglas_ED_Fig5.pdf',
       plot = mobilome_enrichment_by_subtype,
       device = 'pdf',
       dpi = 300,
       width = 6,
       height = 6)


# Write out source data:
source_out <- mobilome_COG_hits[, c('sp', 'partition', 'category', 'Description', 'Subtype',
                                    'OR', 'log2OR', 'p', 'fdr')]

colnames(source_out) <- c('species', 'partition', 'COG_gene_family', 'description', 'mobilome_subtype',
                          'odds_ratio', 'log2_odds_ratio', 'p_value', 'fdr_corrected_p_value')

write.table(x = source_out,
            file = "/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data/ED_Fig5.tsv",
            col.names = TRUE, row.names = FALSE, sep = '\t', quote = FALSE)
