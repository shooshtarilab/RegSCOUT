library(Rsamtools)

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
# converting gene_ids to gene symbols, something along these lines, this is taken directly from hossein's code and 
# has not been modified for use in this Rscript, minor changes may be required

# ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl") 
# 
# gene_ids = c()
# for (i in 1:length(eqtl_results)) {
#   gene_ids = c(gene_ids, eqtl_results[[i]]$gene_id)
# }
# 
# gene_ids = unique(gene_ids)
# attributes_ensemble = listAttributes(ensembl)
# 
# gene_names = getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
#                    filters = "ensembl_gene_id",
#                    values = gene_ids,
#                    mart = ensembl)
# 
# for (i in 1:length(eqtl_results)) {
#   eqtl_results[[i]] = merge(eqtl_results[[i]], gene_names, by.x = "gene_id", by.y = "ensembl_gene_id")
# }
# 
# for (i in 1:length(eqtl_results)) {
#   eqtl_results[[i]] = eqtl_results[[i]][,c("gene_id","external_gene_name", "pvalue", "rsid")]
# }

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





