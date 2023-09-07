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

# Values to report in text:
mean(mobilome_COG_hits[which(mobilome_COG_hits$Subtype == 'Transposon'), 'OR'])
sd(mobilome_COG_hits[which(mobilome_COG_hits$Subtype == 'Transposon'), 'OR'])

mobilome_enrichment_by_subtype <- ggplot(data = mobilome_COG_hits, aes(x = Subtype, y = log2OR)) +
                                        geom_boxplot(outlier.shape = NA) +
                                        geom_quasirandom(alpha = 0.8) +
                                        theme_bw() +
                                        xlab("Associated mobile genetic element sub-type") +
                                        ylab(expression('log'[2]*'(Odd\'s ratio + 0.1)'))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/extended_indepth_mobilome_enrichment.pdf',
       plot = mobilome_enrichment_by_subtype,
       device = 'pdf',
       dpi = 300,
       width = 6,
       height = 6)
