rm(list = ls(all.names = TRUE))

# Data distributions for different taxonomic classes across main dataset.

library(cowplot)
library(ggbeeswarm)
library(ggplot2)
library(ggpubr)

pangenome <- read.table("/data1/gdouglas/projects/pangenome_pseudogene_null_zenodo/broad_pangenome_analysis/pangenome_and_related_metrics_filt.tsv.gz",
                        header = TRUE, sep = "\t", stringsAsFactors = FALSE, row.names = 1)

class_tallies <- table(pangenome$class)
rare_classes <- names(class_tallies)[which(class_tallies <= 5)]
pangenome$class_clean <- pangenome$class
pangenome$class_clean[which(pangenome$class %in% rare_classes)] <- "Other"
pangenome$class_clean <- gsub('c__', '', pangenome$class_clean)
pangenome$class_clean <- factor(pangenome$class_clean,
                                levels = rev(sort(unique(pangenome$class_clean))))

metrics_to_plot <- c('ds', 'dnds',
                     'mean_num_genes', 'mean_num_pseudo',
                     'genomic_fluidity',
                     'mean_percent_singletons_per9', 'mean_percent_singletons_pseudo_per9',
                     'si_sp')

clean_labels <- list('ds' = 'dS',
                     'dnds' = 'dN/dS',
                     'genome_size' = 'Genome size',
                     'mean_num_genes' = 'Mean number of genes',
                     'mean_num_pseudo' = 'Mean number of pseudogenes',
                     'genomic_fluidity' = 'Genomic fluidity',
                     'pseudogene_genomic_fluidity' = 'Genomic pseudo. fluidity',
                     'mean_num_singletons_per9' = 'Mean no. singletons',
                     'mean_num_singletons_pseudo_per9' = 'Mean no. pseudo. singletons',
                     'mean_percent_singletons_per9' = 's[i]',
                     'mean_percent_singletons_pseudo_per9' = 's[p]',
                     'si_sp' = 'si/sp')

all_plots <- list()

metric_count = 1

for (m in metrics_to_plot) {

    all_plots[[m]] <- ggplot(data = pangenome,
                             aes_string(y = 'class_clean',
                                        x = m)) +
                      geom_quasirandom(size = 0.5, orientation = 'y', colour = 'grey70') +
                      geom_boxplot(outlier.shape = NA, fill = 'grey50', alpha = 0.5) +
                      theme_bw() +
                      ylab('')
    
    if (m %in% c('mean_percent_singletons_per9', 'mean_percent_singletons_pseudo_per9', 'si_sp')) {
      all_plots[[m]] <- all_plots[[m]] + xlab(parse(text = clean_labels[m]))
    } else {
      all_plots[[m]] <- all_plots[[m]] + xlab(clean_labels[m])
    }
    
    # Remove y-axis labels for every second panel.
    if (metric_count %% 2 == 0) {
      all_plots[[m]] <- all_plots[[m]] + theme(axis.text.y=element_blank())
    }
    
    metric_count <- metric_count + 1
    
}

# Fix si/sp x-axis label.
all_plots[['si_sp']] <- all_plots[['si_sp']] + xlab(expression(paste(s[i], '/', s[p], sep = '')))

# Put dN/dS plot on log-10 scale and indicate this in x-axis label.
all_plots[['dnds']] <- all_plots[['dnds']] + scale_x_continuous(trans='log10') + xlab(expression('dN/dS (log'[10]*')'))

# Increase x-axis range a bit for dS:
all_plots[['ds']] <- all_plots[['ds']] + xlim(0, 0.22)



combined <- plot_grid(plotlist = all_plots,
                      nrow = 4,
                      rel_widths = c(1, 0.75))

ggsave(filename = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_items/WORKING_metric_distributions_by_taxa_class.pdf',
       plot = combined,
       device = 'pdf',
       dpi = 600,
       width = 9,
       height = 10)
