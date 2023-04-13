rm(list = ls(all.names = TRUE))

COG_description <- read.table("/data1/gdouglas/db/COG_definitions/cog-20.def.tab",
                              row.names = 1, header = FALSE, sep = "\t", quote = "",
                              stringsAsFactors = FALSE, comment.char = "")

ultra.cloud_enrich <- read.table(file = "/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/indepth_10_species_analysis/COG_enrichment_results/ultra.cloud-COG-gene-enrichments.tsv.gz",
                                 header = TRUE, sep = "\t", stringsAsFactors = FALSE)

ultra.cloud_enrich <- ultra.cloud_enrich[which(ultra.cloud_enrich$fdr < 0.05), ]

ultra.cloud_enrich$Descrip <- tolower(COG_description[ultra.cloud_enrich$category, "V3"])


# C - ultra.cloud, pseudogene-enriched, redundant
ultra.cloud_enrich_C_enriched_redundant <- ultra.cloud_enrich[which(ultra.cloud_enrich$COG_category == "C" & ultra.cloud_enrich$type == "pos.enriched_COG_false.non.redundant"), ]

ultra.cloud_enrich_C_enriched_redundant[which(ultra.cloud_enrich_C_enriched_redundant$category == "COG0243"), ]


# F - ultra.cloud, pseudogene-enriched, redundant
ultra.cloud_enrich_F_enriched_redundant <- ultra.cloud_enrich[which(ultra.cloud_enrich$COG_category == "F" & ultra.cloud_enrich$type == "pos.enriched_COG_false.non.redundant"), ]


# J - ultra.cloud, pseudogene-enriched, redundant
ultra.cloud_enrich_J_enriched_redundant <- ultra.cloud_enrich[which(ultra.cloud_enrich$COG_category == "J" & ultra.cloud_enrich$type == "pos.enriched_COG_false.non.redundant"), ]


# Multi-species hits.
ultra.cloud_enrich_J_enriched_redundant[which(ultra.cloud_enrich_J_enriched_redundant$category == "COG0456"), ]
ultra.cloud_enrich_J_enriched_redundant[which(ultra.cloud_enrich_J_enriched_redundant$category == "COG1208"), ]
ultra.cloud_enrich_J_enriched_redundant[which(ultra.cloud_enrich_J_enriched_redundant$category == "COG1418"), ]
ultra.cloud_enrich_J_enriched_redundant[which(ultra.cloud_enrich_J_enriched_redundant$category == "COG1670"), ]

COG_description[c("COG0456", "COG1208", "COG1418", "COG1670"), ]


ultra.cloud_enrich_J_enriched_redundant[which(ultra.cloud_enrich_J_enriched_redundant$OR > 10), ]

# S - ultra.cloud, pseudogene-enriched, redundant
ultra.cloud_enrich_S_enriched_redundant <- ultra.cloud_enrich[which(ultra.cloud_enrich$COG_category == "S" & ultra.cloud_enrich$type == "pos.enriched_COG_false.non.redundant"), ]

multispecies_S_COGs <- names(table(ultra.cloud_enrich_S_enriched_redundant$category))[which(table(ultra.cloud_enrich_S_enriched_redundant$category) > 1)]

ultra.cloud_enrich_S_enriched_redundant[which(ultra.cloud_enrich_S_enriched_redundant$category %in% multispecies_S_COGs), ]


# D - ultra.cloud, pseudogene-depleted, redundant
ultra.cloud_D_depleted_redundant <- ultra.cloud_enrich[which(ultra.cloud_enrich$COG_category == "D" & ultra.cloud_enrich$type == "neg.enriched_COG_false.non.redundant"), ]

ultra.cloud_D_depleted_redundant <- ultra.cloud_D_depleted_redundant[which(ultra.cloud_D_depleted_redundant$OR < 1), ]

multispecies_D_COGs <- names(table(ultra.cloud_D_depleted_redundant$category))[which(table(ultra.cloud_D_depleted_redundant$category) > 1)]

COG_description[multispecies_D_COGs, ]

ultra.cloud_D_depleted_redundant[which(ultra.cloud_D_depleted_redundant$category %in% multispecies_D_COGs), ]

ultra.cloud_D_depleted_redundant[which(ultra.cloud_D_depleted_redundant$category == "COG1192"), ]