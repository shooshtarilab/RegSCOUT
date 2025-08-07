library(Rsamtools)
library(biomaRt)
library(dplyr)

# read in tabix dataset
eqtl_dir <- '/Users/richardzhang/Desktop/OneDrive/Desktop/Thesis_Project/RegSCOUT_github_update/eqtl_test_data/QTD000115.cc.tsv.gz'

eqtl_df <- read.table(
  gzfile(eqtl_dir),
  header = TRUE,
  sep = "\t"
)

# make necessary modifications
eqtl_df <- eqtl_df[eqtl_df$pvalue < 0.001,]
eqtl_df$chromosome <- paste0('chr', eqtl_df$chromosome)
eqtl_df <- eqtl_df[,c('chromosome', 'position', 'gene_id', 'rsid')]

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
  file = "/Users/richardzhang/Desktop/OneDrive/Desktop/Thesis_Project/RegSCOUT_github_update/eqtl_test_data/test_data.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

# adding the column names back in, but as a commented line (i.e., starting with #)
writeLines(
  c(column_names, readLines("/Users/richardzhang/Desktop/OneDrive/Desktop/Thesis_Project/RegSCOUT_github_update/eqtl_test_data/test_data.tsv")),
  con = "/Users/richardzhang/Desktop/OneDrive/Desktop/Thesis_Project/RegSCOUT_github_update/eqtl_test_data/test_data.tsv"
)

# run these lines of code in command line make sure bgzip and tabix are installed
# bgzip < test_data.tsv > test_data.tsv.gz (first file is input, second file is output)
# tabix -s 1 -b 2 -e 2 test_data.tsv.gz
  # -s is the number of the column that has chromosome numbers
  # -b is the column that has start position
  # -e is the end position
  # -b and -e can be identical





