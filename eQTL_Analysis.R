library(dplyr)
library(stringr)
library(R.utils)
library(GenomicRanges)
library(Rsamtools)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)
# read output directory
output_dir = args[["output_dir"]]

# read in SNP, rmp information
snp_rmp_dir <- paste0(output_dir, "risk_regions_ratio.txt") # rmp region, cell type adn TFSNP
snp_rmp_df = read.table(paste0(output_dir, "risk_regions_ratio.txt"), header = TRUE)
# modify this dataframe
snp_rmp_df$log_lik_ratio <- NULL
colnames(snp_rmp_df)[colnames(snp_rmp_df) == "region"] <- "rmp"

# separating snps from TFs in TFSNP column
snp_rmp_df$snp <- gsub(".*-", "", snp_rmp_df$TFSNP)
snp_rmp_df$TFSNP <- NULL
snp_rmp_df <- snp_rmp_df %>% distinct()

# adding snp position and chromosome information
effect_snp_results <- read.table(paste0(output_dir, "Ci_effect_SNPs.txt"), header = T, sep = '\t')
effect_snp_results <- effect_snp_results[,c('SNP', 'CHR', 'Pos')] %>% distinct()

snp_rmp_df <- snp_rmp_df %>%
  left_join(
    effect_snp_results,
    by = c("snp" = "SNP")
  )

genome_build = args[["genome_built"]]

# read in user instructions
eqtl_instruct_dir <- args[["eqtl_instruct_dir"]]
user_instruct <- read.table(eqtl_instruct_dir, header = TRUE)

# creating list of results for eQTL
eqtl_results_list <- list()

# iterating over number of rows of user_instruct
num_eqtl <- nrow(user_instruct)

count = 1 # needed
for (i in 1:num_eqtl) {
  current_row <- user_instruct[i,]
  eqtl_dir <- current_row$eqtl_dir
  tabix <- current_row$tabix
  atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
  
  if (!as.logical(tabix)) {
    # read in eqtl dataset
    eqtl_dataset <- read.table(file = eqtl_dir, sep = '\t', header = TRUE)
    
    # filter snp_rmp_df
    snp_rmp_filt <- snp_rmp_df[snp_rmp_df$cell %in% atac_cell_types,]
    
    if ('chr' %in% colnames(eqtl_dataset) & 'pos' %in% colnames(eqtl_dataset)) { # if chromosome and position columns are present then using those columns
      # preparing for inner_join between snp_rmp_filt and eqtl_dataset
      snp_rmp_filt$Pos <- as.integer(snp_rmp_filt$Pos)
      eqtl_dataset$pos <- as.integer(eqtl_dataset$pos)
      eqtl_dataset$snp <- NULL
      
      # merge eqtl with snp_rmp dataframe
      eqtl_results <- snp_rmp_filt %>%
        inner_join(
          eqtl_dataset,
          by = c("CHR" = "chr", "Pos" = "pos")
        )
    } else { # using snp column if chr and pos columns are not provided
      # merge eqtl with snp_rmp dataframe
      eqtl_results = merge(snp_rmp_filt, eqtl_dataset, by = 'snp')
    }
    
    # clean up dataframe and append to list
    eqtl_results <- eqtl_results[,c('snp', 'rmp', 'cell', 'gene')]
    eqtl_results_list[[i]] <- eqtl_results
    
    # save the dataframe
    write.table(eqtl_results, file = paste0(output_dir, "eqtl_analysis", i, "_results.txt"), row.names = F, quote = F,
                sep = '\t')
  } else { # if tabix files are provided
    # read in eqtl dataset
    eqtl_dataset <- TabixFile(eqtl_dir)
    
    for (cell_type in atac_cell_types){
      # filter snp_rmp_df
      snp_rmp_filt <- snp_rmp_df[snp_rmp_df$cell %in% cell_type,]

      # defining genomic ranges for SNPs
      snp_granges <- GRanges(seqnames = snp_rmp_filt$CHR,
                                ranges = IRanges(start = snp_rmp_filt$Pos, end = snp_rmp_filt$Pos))
      
      # snp_granges@seqnames@values <- gsub("^chr", "", snp_granges@seqnames@values)
      if (genome_build == "hg19"){
        hg19to38 = import.chain(paste0(output_dir, "hg19ToHg38.over.chain"))
        snp_granges = unlist(liftOver(snp_granges,hg19to38))
      }
      eqtl_tabix <- scanTabix(eqtl_dataset, param = snp_granges)
      
      # add metadata from snp_rmp_filt for each matched row
      eqtl_results_meta <- lapply(seq_along(eqtl_tabix), function(j) {
        lines <- eqtl_tabix[[j]]
        split_lines <- strsplit(lines, "\t")
        eqtl_results_df <- as.data.frame(do.call(rbind, split_lines), stringsAsFactors = FALSE)
        meta_row <- snp_rmp_filt[j, , drop = FALSE]
        add_meta <- meta_row[rep(1, nrow(eqtl_results_df)), , drop = FALSE]
        cbind(add_meta, eqtl_results_df)
      })
      
      # create a dataframe and save it in list
      eqtl_results_df <- do.call(rbind, eqtl_results_meta)
      eqtl_results_df <- eqtl_results_df[,c('snp', 'rmp', 'cell', 'V3')] 

      colnames(eqtl_results_df)[colnames(eqtl_results_df) == 'V3'] <- 'gene'
      eqtl_results_list[[count]] <- eqtl_results_df
      # save the dataframe
      write.table(eqtl_results_df, file = paste0(output_dir, "eqtl_analysis", count,"_",cell_type, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
      count = count + 1
    }
  }
}

# obtaining one dataframe with all results
all_results <- bind_rows(eqtl_results_list)  
rownames(all_results) <- NULL

# save this dataframe
write.table(all_results, file = paste0(output_dir, "all_eqtl_results.txt"), row.names = F, quote = F,
            sep = '\t')

print('eQTL analysis complete!')
  