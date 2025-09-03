suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(biomaRt))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(data.table))

# read in tabix dataset
eqtl_dir <- '/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/eqtl_files'

eqtl_files <- list.files(eqtl_dir, pattern = "\\.gz$", full.names = TRUE)

# Iterate over them
for (file_path in eqtl_files) {
  file_path = "/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/eqtl_files/QTD000066.all.tsv.gz"
  filename = basename(file_path)
  filename <- sub("\\.gz$", "", filename)
  print(filename)
  # Read only necessary columns
  eqtl_df <- fread(
    file_path,
    sep = "\t",
    header = TRUE,
    select = c("chromosome", "position", "gene_id", "rsid", "pvalue")
  )

  # Filter and transform
  eqtl_df <- eqtl_df[pvalue < 0.001]
  eqtl_df[, chromosome := paste0("chr", chromosome)]
  eqtl_df <- eqtl_df[, .(chromosome, position, gene_id, rsid)]

  # once the dataframe has been filtered for significance and certain columns, biomaRt can be used to bring in gene names 
  # from ensembl ids if desired
  ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") # sometimes ensembl servers are busy and this function won't work, usually waiting a bit before running this function again helps there

  # get list of ensembl gene ids from eqtl dataframe
  gene_ids = unique(c(eqtl_df$gene_id))

  # getting gene names using biomaRt
  gene_names = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id",
                    values = gene_ids,
                    mart = ensembl)
  gene_names$external_gene_name[gene_names$external_gene_name == ''] <- NA # replace blank gene names with NA

  # adding gene names to eqtl_df, for those ensembl ids where gene names not found, keeping the ensembl id
  eqtl_df <- eqtl_df %>%
    left_join(
      gene_names,
      by = c("gene_id" = "ensembl_gene_id")
    )

  eqtl_df$external_gene_name[is.na(eqtl_df$external_gene_name)] <- eqtl_df$gene_id[is.na(eqtl_df$external_gene_name)]
  eqtl_df <- eqtl_df[,c('chromosome', 'position', 'external_gene_name', 'rsid')] %>% distinct()

  # storing column names in variable
  column_names <- paste("#", paste(colnames(eqtl_df), collapse = "\t"))

  # save this as a tsv file
  write.table(
    eqtl_df,
    file = paste0("/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/eqtl_filtered_files/",filename),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )

  # adding the column names back in, but as a commented line (i.e., starting with #)
  writeLines(
    c(column_names, readLines(paste0("/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/eqtl_filtered_files/",filename))),
    con = paste0("/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/eqtl_filtered_files/",filename)
  )
  input_file = paste0("/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/eqtl_filtered_files/",filename)
  output_file <- paste0(input_file,".gz")

  # Compress with bgzip
  system(paste("bgzip -f", shQuote(input_file)))

  # Index with tabix
  system(paste("tabix -s 1 -b 2 -e 2", shQuote(output_file)))
  print(c("done ", filename))
}



# run these lines of code in command line make sure bgzip and tabix are installed
# bgzip < test_data.tsv > test_data.tsv.gz (first file is input, second file is output)
# tabix -s 1 -b 2 -e 2 test_data.tsv.gz
  # -s is the number of the column that has chromosome numbers
  # -b is the column that has start position
  # -e is the end position
  # -b and -e can be identical





