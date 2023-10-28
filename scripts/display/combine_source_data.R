rm(list = ls(all.names = TRUE))

# Combine source data-files into spreadsheet.
library(openxlsx)

source_data <- list()

# Read in all source data-files to workbook.
main_source <- list.files(path = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data',
                          pattern = '^Fig.*.tsv.gz',
                          full.names = TRUE)

extended_source <- list.files(path = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data',
                              pattern = '^ED_.*.tsv.gz',
                              full.names = TRUE)

all_source <- c(main_source, extended_source)

combined_source_data <- createWorkbook()

sheet_num <- 1

for (source_file in all_source) {

  source_df <- read.table(file = source_file,
                          header = TRUE, sep = '\t', stringsAsFactors = FALSE)

  source_id <- gsub('.tsv.gz$', '', basename(source_file))

  addWorksheet(combined_source_data, sheetName = source_id)

  writeData(combined_source_data, sheet = sheet_num, x = source_df)
  
  sheet_num <- sheet_num + 1
}

saveWorkbook(wb = combined_source_data,
             file = '/home/gdouglas/scripts/pangenome_pseudogene_null/display_source_data/Douglas_display_source_data.xlsx')

# Note that some column headers were then manually cleaned up.
# ED Table 1 was especially cleaned up to remove extra information.
