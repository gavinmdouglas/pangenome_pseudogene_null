rm(list = ls(all.names = TRUE))

### Parse eggNOG annotation and define COG categories based on the presence of *any* COG amongst list of OGs that were hit.

eggNOG_annot <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/cdhit_workflow/eggnog_mapper_annot/derep_gene_annot.emapper.annotations.tsv.gz",
                           header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = FALSE, quote = "", comment.char = "")

eggNOG_annot$all_COG <- sapply(eggNOG_annot$eggNOG_OGs, function(x) { paste(grep("^COG", unique(gsub("@.*$", "", strsplit(x, ",")[[1]])), value = TRUE), collapse = ",") })
  
COG2catgory_uniq <- readRDS("/data1/gdouglas/db/COG_definitions/cog-20.to_category_collapse.rds")
  
rep2cluster <- read.table("/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/cdhit_workflow/cluster2rep.tsv.gz",
                            header = FALSE, sep = "\t", row.names = 2, stringsAsFactors = FALSE)
  
COG_annot <- eggNOG_annot[, c("all_COG", "Description"), drop = FALSE]

rownames(COG_annot) <- rep2cluster[rownames(COG_annot), "V1"]

COG_annot$COG_category <- "-"

clusters_w_any.COG <- grep("^COG", COG_annot$all_COG)

COG_annot$COG_category[clusters_w_any.COG] <- sapply(clusters_w_any.COG,
                                                              function(i) {
                                                    
                                                                all_COGs <- COG_annot[i, "all_COG"]
                                                                categories <- c()
                                                                COGs <- strsplit(all_COGs, ",")[[1]]
                                                                for (COG in COGs) {
                                                                  categories <- c(categories, COG2catgory_uniq[COG, "category"])
                                                                }
                                                                
                                                                return(paste(sort(unique(categories)), collapse = ","))
                                                              })

COG_annot$COG_category[which(COG_annot$COG_category == "")] <- "-"

saveRDS(object = COG_annot,
        file = "/data1/gdouglas/projects/hgt_fragments/pseudogenes/processed/cdhit_workflow/eggnog_mapper_annot/cluster_annot.rds")
