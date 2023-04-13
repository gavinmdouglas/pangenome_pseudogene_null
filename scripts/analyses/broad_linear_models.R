rm(list = ls(all.names = TRUE))

# Fit linear models on pangenome diversity metrics vs molecular evolution metrics and taxonomic classes.

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

# Classes with < 5 members were collapsed into "Other". These genomes were used as the first level of the class variable, i.e., used for the intercept.
class_tallies <- table(pangenome$class)
rare_classes <- names(class_tallies)[which(class_tallies < 5)]
pangenome$class_clean <- pangenome$class
pangenome$class_clean[which(pangenome$class %in% rare_classes)] <- "Other"
class_order <- sort(unique(pangenome$class_clean))
pangenome$class_clean <- factor(pangenome$class_clean, levels = c("Other", class_order[-which(class_order == "Other")]))

mean.percent.ratio_working <- data.frame(si_sp = log10(pangenome$si_sp),
                                         log_c_dS = log10(pangenome$median_pairwise_ds),
                                         log_c_dNdS = log10(pangenome$median_MG94_dnds),
                                         class = pangenome$class_clean)

# Mean-centre log dS and dNdS values.
mean.percent.ratio_working$log_c_dS <- mean.percent.ratio_working$log_c_dS - mean(mean.percent.ratio_working$log_c_dS)
mean.percent.ratio_working$log_c_dNdS <- mean.percent.ratio_working$log_c_dNdS - mean(mean.percent.ratio_working$log_c_dNdS)

mean.percent.ratio_models <- list()
mean.percent.ratio_models[['class']] <- lm(si_sp ~ class, data = mean.percent.ratio_working)
mean.percent.ratio_models[['class_dS']] <- lm(si_sp ~ class + log_c_dS, data = mean.percent.ratio_working)
mean.percent.ratio_models[['class_dNdS']] <- lm(si_sp ~ class + log_c_dNdS, data = mean.percent.ratio_working)
mean.percent.ratio_models[['class_dS_dNdS']] <- lm(si_sp ~ class + log_c_dS + log_c_dNdS, data = mean.percent.ratio_working)



mean.percent_working <- data.frame(singletons_percent = log10(pangenome$mean_percent_singletons_per9),
                                   log_c_dS = log10(pangenome$median_pairwise_ds),
                                   log_c_dNdS = log10(pangenome$median_MG94_dnds),
                                   class = pangenome$class_clean)

# Mean-centre log dS and dNdS values.
mean.percent_working$log_c_dS <- mean.percent_working$log_c_dS - mean(mean.percent_working$log_c_dS)
mean.percent_working$log_c_dNdS <- mean.percent_working$log_c_dNdS - mean(mean.percent_working$log_c_dNdS)

mean.percent_models <- list()
mean.percent_models[['class']] <- lm(singletons_percent ~ class, data = mean.percent_working)

mean.percent_models[['class_dS']] <- lm(singletons_percent ~ class + log_c_dS, data = mean.percent_working)

mean.percent_models[['class_dNdS']] <- lm(singletons_percent ~ class + log_c_dNdS, data = mean.percent_working)

mean.percent_models[['class_dS_dNdS']] <- lm(singletons_percent ~ class + log_c_dS + log_c_dNdS, data = mean.percent_working)

out_models <- list(si_sp = mean.percent.ratio_models,
                   si = mean.percent_models)

saveRDS(object = out_models,
        file = '/data1/gdouglas/projects/pangenome_pseudogene_null_figshare/broad_pangenome_analysis/model_output/pangenome_linear_models.rds')
