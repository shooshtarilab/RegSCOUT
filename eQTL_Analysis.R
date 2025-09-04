suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Rsamtools))
suppressPackageStartupMessages(library(liftOver))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)
# read output directory
output_dir = args[["output_dir"]]

genome_build = args[["genome_build"]]

# read in SNP, rmp information
snp_rmp_dir <- paste0(output_dir, "risk_regions_ratio.txt")
snp_rmp_df = read.table(paste0(output_dir, "risk_regions_ratio.txt"), header = TRUE)
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

# read in user instructions
eqtl_instruct_dir <- args[["eqtl_instruct_dir"]]
if (!file.exists(eqtl_instruct_dir)) {
  stop("eQTL analysis requested but eQTL instructions spreadsheet not found, please ensure path is correct/provided. Or if eQTL analysis is not desired please do not use --hic_eqtl_analysis parameter.")
} 
user_instruct <- read.table(eqtl_instruct_dir, header = TRUE)

# getting list of atac cell types requested and rmp cell types, seeing if rmps were not found in some atac cell types, removing them
rmp_cell_types <- unique(snp_rmp_df$cell)
peak_cell_types <- unique(unlist(str_split(user_instruct$atac_cell_types, ',')))
missing_cell_types <- setdiff(peak_cell_types, rmp_cell_types)

if (length(missing_cell_types) > 0) { # removing missing cell types from user instructions
  user_instruct <- user_instruct %>%
    rowwise() %>%
    mutate(
      atac_cell_types = str_split(atac_cell_types, ",")[[1]] %>%
        setdiff(missing_cell_types) %>%
        paste(collapse = ",")
    ) %>%
    filter(atac_cell_types != "") %>%
    ungroup()
  
  message("Note: eQTL analysis will not be conducted on cell types in which no SNPs co-localizing with risk-mediating peaks were identified: ", paste(missing_cell_types, collapse = ", "), ". Analysis not being conducted may also occur if a cell type name in the atac_cell_types column does not match with the corresponding cell type name provided in the scATAC-seq object (i.e., B cell vs b_cell).")
}

# check to see if any instructions entries are left
if (nrow(user_instruct) == 0) {
  stop("No eQTL instruction entries remain after removing ATAC cell types requested but not found in risk-mediating peak results. eQTL analysis will not be conducted.")
}

# creating list of results for eQTL
eqtl_results_list <- list()

# iterating over number of rows of user_instruct
num_eqtl <- nrow(user_instruct)

for (i in 1:num_eqtl) {
  current_row <- user_instruct[i,]
  eqtl_dir <- current_row$eqtl_dir
  tabix <- as.logical(current_row$tabix)
  atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
  
  if (!tabix) {
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
    
    # see if any results were found
    if (nrow(eqtl_results) == 0) {
      message("Note: No eQTL results found for row with eQTL directory '", eqtl_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. One possible reason for this is a mismatch between genome builds of eQTL and SNP data.")
    } else {
      # clean up dataframe and append to list
      eqtl_results <- eqtl_results[,c('snp', 'rmp', 'cell', 'gene')]
      eqtl_results_list[[i]] <- eqtl_results
      
      # save the dataframe
      write.table(eqtl_results, file = paste0(output_dir, "eqtl_analysis", i, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
    }
  } else if (tabix) { # if tabix files are provided
    # read in eqtl dataset
    eqtl_dataset <- TabixFile(eqtl_dir)
    
    # filter snp_rmp_df
    snp_rmp_filt <- snp_rmp_df[snp_rmp_df$cell %in% atac_cell_types,]
    

    # defining genomic ranges for SNPs
    snp_granges <- GRanges(seqnames = snp_rmp_filt$CHR,
                               ranges = IRanges(start = snp_rmp_filt$Pos, end = snp_rmp_filt$Pos))
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
    
    if (nrow(eqtl_results_df) == 0) {
      message("Note: No eQTL results found for row with eQTL directory '", eqtl_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. One possible reason for this is a mismatch between genome builds of eQTL and SNP data.")
    } else {
      eqtl_results_df <- eqtl_results_df[,c('snp', 'rmp', 'cell', 'V3')] 

      colnames(eqtl_results_df)[colnames(eqtl_results_df) == 'V3'] <- 'gene'
      eqtl_results_list[[i]] <- eqtl_results_df
      
      # save the dataframe
      write.table(eqtl_results_df, file = paste0(output_dir, "eqtl_analysis", i, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
    }
  }
}

# check to see if any results were obtained
if (length(eqtl_results_list) == 0) {
  message('eQTL analysis complete, no results were found over all eQTL datasets.')
} else {
  # obtaining one dataframe with all results
  all_results <- bind_rows(eqtl_results_list) %>% distinct()
  rownames(all_results) <- NULL
  
  # save this dataframe
  write.table(all_results, file = paste0(output_dir, "all_eqtl_results.txt"), row.names = F, quote = F,
              sep = '\t')
  
  message('eQTL analysis complete!')
}

  