library(readxl)
library(dplyr)
library(stringr)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# read output directory
output_dir = args[["output_dir"]]

# read in SNP, rmp information
snp_rmp_dir <- paste0(output_dir, "risk_regions_ratio.xlsx")
snp_rmp_df <- read_xlsx(snp_rmp_dir)

# modify this dataframe
snp_rmp_df$log_lik_ratio <- NULL
colnames(snp_rmp_df)[colnames(snp_rmp_df) == "region"] <- "rmp"

# separating snps from TFs in TFSNP column
snp_rmp_df$snp <- gsub(".*-", "", snp_rmp_df$TFSNP)
snp_rmp_df$TFSNP <- NULL
snp_rmp_df <- snp_rmp_df %>% distinct()

# read in user instructions
eqtl_instruct_dir <- args[["eqtl_instruct_dir"]]
user_instruct <- read_xlsx(eqtl_instruct_dir)

# creating list of results for Hi-C
eqtl_results_list <- list()

# iterating over number of rows of user_instruct
num_eqtl <- nrow(user_instruct)

for (i in 1:num_eqtl) {
  current_row <- user_instruct[i,]
  eqtl_dir <- current_row$eqtl_dir
  
  # read in eqtl dataset and associated information
  eqtl_dataset <- read.table(file = eqtl_dir, sep = '\t', header = TRUE)
  atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
  
  # filter snp_rmp_df
  snp_rmp_filt <- snp_rmp_df[snp_rmp_df$cell %in% atac_cell_types,]
  
  # merge eqtl with snp_rmp dataframe
  eqtl_results = merge(snp_rmp_filt, eqtl_dataset, by = 'snp')
  
  # clean up dataframe and append to list
  eqtl_results <- eqtl_results[,c('snp', 'rmp', 'cell', 'gene')]
  eqtl_results_list[[i]] <- eqtl_results
  
  # save the dataframe
  write.table(eqtl_results, file = paste0(output_dir, "eqtl_analysis", i, "_results.txt"), row.names = F, quote = F,
              sep = '\t')
}

# obtaining one dataframe with all results
all_results <- bind_rows(eqtl_results_list)  
rownames(all_results) <- NULL

# save this dataframe
write.table(all_results, file = paste0(output_dir, "all_eqtl_results.txt"), row.names = F, quote = F,
            sep = '\t')

print('eQTL analysis complete!')
  