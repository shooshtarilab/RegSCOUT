suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# read output directory
output_dir = args[["output_dir"]]

# defining default parameter values
defaults <- list(
  prom_th_up = 2000,
  prom_th_down = 2000,
  quantile_th = 0.25
)

# read in TF, SNP, and rmp information
tf_table <- read.table(paste0(output_dir, "risk_regions_ratio.txt"), header = TRUE)
tf_table <- tf_table %>%
  mutate(
    tf = str_extract(TFSNP, "^.+(?=-[^-]+$)"),
    snp = str_extract(TFSNP, "[^\\-]+$") 
  )

# just keep TF and cell columns
tf_table_filt <- tf_table %>%
  select(tf, cell) %>%
  distinct()

# defining function that identifies peaks on TF promoters
confirm_tf_promoter_peaks <- function(tf_list, heterodimer_list, prom_th_up, prom_th_down, peak_file_path, gencode_file_path, prio_tf_table) {
  gene_tss_grg = readRDS(paste0(output_dir, "gene_tss_granges.rds"))

  gene_tss_grg = gene_tss_grg[gene_tss_grg$gene_name %in% tf_list, ]

  # Read peak data
  peaks <- read.delim(peak_file_path, header = TRUE)
  peaks_ranges <- IRanges(start = as.integer(peaks$start), end = as.integer(peaks$end))
  peaks_granges <- GRanges(seqnames = peaks$chr, ranges = peaks_ranges, cell = peaks$cell_sub_types)
  
  # Find overlaps
  TF_peak_overlap <- findOverlaps(gene_tss_grg, peaks_granges)
  
  if (length(TF_peak_overlap) == 0) {
    return('none')
  } else {
    TF_peak_overlap_df <- as.data.frame(TF_peak_overlap)
    TF_peak_overlap_df$gene_name <- gene_tss_grg$gene_name[TF_peak_overlap_df$queryHits]
    TF_peak_overlap_df$gene_chr <- as.character(seqnames(gene_tss_grg)[TF_peak_overlap_df$queryHits])
    TF_peak_overlap_df$gene_start <- start(gene_tss_grg)[TF_peak_overlap_df$queryHits]
    TF_peak_overlap_df$gene_end <- end(gene_tss_grg)[TF_peak_overlap_df$queryHits]
    TF_peak_overlap_df$peak_chr <- as.character(seqnames(peaks_granges)[TF_peak_overlap_df$subjectHits])
    TF_peak_overlap_df$peak_start <- start(peaks_granges)[TF_peak_overlap_df$subjectHits]
    TF_peak_overlap_df$peak_end <- end(peaks_granges)[TF_peak_overlap_df$subjectHits]
    TF_peak_overlap_df$peak_cell <- peaks_granges$cell[TF_peak_overlap_df$subjectHits]
    
    # simplify overlap dataframe
    tf_peak <- TF_peak_overlap_df %>%
      select(gene_name, peak_cell) %>%
      separate_rows(peak_cell, sep = ",") %>%
      distinct()
    colnames(tf_peak) <- c('tf', 'cell')
    
    tf_peak_complete <- tf_peak
    
    # account for TF heterodimers
    if (length(heterodimer_list) > 0) {
      heterodimer_parts <- strsplit(heterodimer_list, "::")
      
      for (i in 1:length(heterodimer_parts)) {
        current_hd <- heterodimer_parts[[i]]
        tf1 <- current_hd[1]
        tf2 <- current_hd[2]
        tf_peak_df1 <- tf_peak[tf_peak$tf == tf1,]
        tf_peak_df2 <- tf_peak[tf_peak$tf == tf2,]
        if (nrow(tf_peak_df1) > 0 & nrow(tf_peak_df2) > 0) {
          cell_overlap <- intersect(tf_peak_df1$cell, tf_peak_df2$cell)
          hd_row <- as.data.frame(matrix(nrow = length(cell_overlap), ncol = 2))
          colnames(hd_row) <- c('tf', 'cell')
          hd_row$tf <- heterodimer_list[i]
          hd_row$cell <- cell_overlap
          # adding this information to tf_peak_complete
          tf_peak_complete <- rbind(tf_peak_complete, hd_row)
        }
      }
    }
    
    # now identifying which cell type-specific TFs identified previously are also potentially expressed in the cell type
    expr_prio_tf <- dplyr::intersect(prio_tf_table, tf_peak_complete)
    
    # check if the result of dplyr intersect is an empty data frame
    if (nrow(expr_prio_tf) == 0) {
      return('none')
    } else {
      return(expr_prio_tf)
    }
  }
}

# defining multiple functions to conduct scRNA-seq analysis for TF expression
# function 1a: calculating percent expression of each gene in each cell type, this is used when scRNA-seq seurat
# has only one cell type
gene_percent_expr_ct <- function(seurat_obj, matrix_location) { 
  # obtaining gene expression data
  assay = seurat_obj[["RNA"]]
  expr_data = slot(assay, matrix_location)
  
  # calculating gene expression frequencies
  results <- rowSums(expr_data > 0) / ncol(expr_data)
  
  results_df <- data.frame(results)
  
  # replace all NAs with 0
  results_df[is.na(results_df)] <- 0
  
  return(results_df)
}

# function 1b: creating a function to create dataframe that lists all genes and the percent of cells in which 
# they are expressed, for each cell type
gene_percent_expr <- function(seurat_obj, matrix_location) {
  cell_types <- levels(seurat_obj@active.ident)
  results_df <- as.data.frame(matrix(nrow = nrow(seurat_obj), ncol = length(cell_types)))
  row.names(results_df) <- row.names(seurat_obj)
  colnames(results_df) <- cell_types
  
  i = 1
  for (cell in cell_types) {
    cell_labels_filt = WhichCells(object = seurat_obj, idents = cell)
    expr_data = FetchData(seurat_obj, row.names(seurat_obj), cells = cell_labels_filt, layer = matrix_location)
    
    # for loop to iterate through each gene
    for (gene in colnames(expr_data)) {
      temp_list = expr_data[,gene]
      times_expr = sum(temp_list > 0)
      results_df[gene,i] = times_expr/nrow(expr_data)
    }
    
    i = i + 1
  }
  
  # replace all NAs with 0
  results_df[is.na(results_df)] <- 0
  
  return(results_df)
}

# function 2: determining the 25% quantile of the fraction of cells that express a gene in a cell type
quant_per_celltype <- function(percent_df, quantile_threshold) {
  cell_types <- colnames(percent_df)
  
  # initiate a vector for 25% quantile values
  quant_vec <- c()
  
  for (cell in cell_types) {
    cell_data <- percent_df[[cell]]
    cell_data <- cell_data[cell_data != 0] # remove all zero values
    quant <- quantile(cell_data, probs = quantile_threshold)
    quant_vec <- c(quant_vec, quant)
  }
  
  names(quant_vec) <- cell_types
  
  return(quant_vec)
}

# function 3a: creating function for TF expression analysis, if scRNA-seq dataset only has one cell type
tf_expression_analysis_ct <- function(seurat_obj, TF_list, quant_val, cell_types, matrix_location) {
  # obtaining gene expression data
  assay = seurat_obj[["RNA"]]
  expr_data = slot(assay, matrix_location)
  
  # obtain TF expression data
  expr_data = expr_data[rownames(expr_data) %in% TF_list,]
  
  # check to see if expression data has not been filtered out due to TF_list
  if (nrow(expr_data) == 0) {
    return('none')
  } else {
    # calculating gene expression frequencies
    results <- rowSums(expr_data > 0) / ncol(expr_data)
    
    # changing values to either 0 (if % expression does not pass quantile threshold) or 1 (if it does)
    threshold = quant_val
    results_thr <- as.numeric(results > threshold)
    names(results_thr) <- names(results)
    
    # check to see that at least one TF passes threshold
    if (sum(results_thr) == 0) {
      return('none')
    } else {
      # creating dataframe
      num_col <- length(cell_types)
      row_names <- names(results)
      results_df <- data.frame(matrix(rep(results_thr, times = num_col), ncol = num_col))
      colnames(results_df) <- cell_types
      rownames(results_df) <- row_names
      results_df <- results_df[rowSums(results_df) > 0, , drop = F] 
      
      results_df <- rownames_to_column(results_df, var = "TFs_prio")
      
      return(results_df)
    }
  }
}

# function 3b: creating function for TF expression analysis, use function based on structure of seurat
tf_expression_analysis <- function(seurat_obj, TF_list, quant_vector, matrix_location) { 
  cell_types <- levels(seurat_obj@active.ident)
  results_df <- as.data.frame(matrix(nrow = length(TF_list), ncol = length(cell_types)))
  row.names(results_df) <- TF_list
  colnames(results_df) <- cell_types
  
  # ensure TFs in TF list are present in seurat object
  tf_overlap = intersect(rownames(seurat_obj), TF_list)
  
  if (length(tf_overlap) == 0) {
    return('none')
  } else {
    # for loop to fill results_df
    i = 1
    for (cell in cell_types) {
      cell_labels_filt = WhichCells(object = seurat_obj, idents = cell)
      
      expr_data = FetchData(seurat_obj, TF_list, cells = cell_labels_filt, layer = matrix_location)
      
      # for loop to iterate through each gene
      for (gene in colnames(expr_data)) {
        temp_list = expr_data[,gene]
        times_expr = sum(temp_list > 0)
        threshold = quant_vector[cell] # defining a threshold based on desired quantile
        if (times_expr/nrow(expr_data) > threshold){ 
          results_df[gene,i] = 1
        } else {
          results_df[gene,i] = 0
        }
      }
      
      i = i + 1
    }
    
    # remove rows where all NA in results dataframe (only occurs if TF not found), and all TFs that were not expressed in any cell type
    filt_results_df <- results_df[!rowSums(is.na(results_df)) == ncol(results_df),]
    filt_results_df <- filt_results_df[rowSums(filt_results_df) > 0,] 
    
    if (nrow(filt_results_df) == 0) {
      return('none')
    } else {
      filt_results_df <- rownames_to_column(filt_results_df, var = "TFs_prio")
      return(filt_results_df)
    }
  }
}

# function 4: creating a function to filter these TF heterodimers such that only those that are expressed are left
# and adding TF heterodimers to expression dataframe
heterodimers_expr <- function(heterodimer_list, expr_df) {
  expr_TFs <- expr_df$TFs_prio
  heterodimer_parts <- strsplit(heterodimer_list, "::")
  
  # for a TF heterodimer to be considered expressed both its protein subunits should be expressed
  expressed_hd <- c()
  for (i in 1:length(heterodimer_parts)) {
    current_hd <- heterodimer_parts[[i]]
    if (all(current_hd %in% expr_TFs)) {
      expressed_hd <- c(expressed_hd, heterodimer_list[i])
    }
  }
  
  # Creating a dataframe for expressed heterodimers and determining cell types they are expressed in
  hd_expr_df <- as.data.frame(matrix(nrow = 0, ncol = ncol(expr_df)))
  colnames(hd_expr_df) <- colnames(expr_df)
  
  for (hd in expressed_hd) {
    TFs_in_hd <- strsplit(hd, "::")[[1]]
    temp_df <- expr_df[expr_df$TFs_prio %in% TFs_in_hd,]
    hd_df <- temp_df %>%
      summarise(
        across(
        where(is.numeric), 
        ~ if(sum(.) == 2) 1 else 0
        )
      )
    hd_df <- hd_df %>%
      mutate(TFs_prio = hd) %>%
      select(TFs_prio, everything())
    hd_expr_df <- rbind(hd_expr_df, hd_df)
  }
  
  # add this heterodimer dataframe to original expression dataframe
  new_df <- rbind(expr_df, hd_expr_df)
  new_df <- column_to_rownames(new_df, var = "TFs_prio")
  new_df <- new_df[rowSums(new_df) > 0, , drop = F] 
  
  return(new_df)
}

# function 5a: identify TFs that are both expressed in a cell type and prioritized in a cell type previously, this function
# is specific to when the scRNA-seq dataset only encompassed one cell type
finalize_tfs_prio_expr_ct <- function(prio_tf_table, expr_tf_table, atac_ct, rna_ct) {
  # turn expr_tf_table into same format as prio_tf_table ('tf' and 'cell' columns)
  expr_tf_table <- rownames_to_column(expr_tf_table, var = "tf")
  expr_tf_table_new <- expr_tf_table %>%
    pivot_longer(
      cols = -tf,
      names_to = "cell",
      values_to = "value"
    ) %>%
    filter(value == 1) %>%
    select(tf, cell)
  
  # now identifying which cell type-specific TFs identified previously are also being transcribed in the cell type
  expr_prio_tf <- dplyr::intersect(prio_tf_table, expr_tf_table_new)
  
  return(expr_prio_tf)
}

# function 5b: identify TFs that are both expressed in a cell type and prioritized in a cell type previously
finalize_tfs_prio_expr <- function(prio_tf_table, expr_tf_table, atac_ct, rna_ct) {
  # turn expr_tf_table into same format as prio_tf_table ('tf' and 'cell' columns)
  expr_tf_table <- rownames_to_column(expr_tf_table, var = "tf")
  expr_tf_table_new <- expr_tf_table %>%
    pivot_longer(
      cols = -tf,
      names_to = "cell",
      values_to = "value"
    ) %>%
    filter(value == 1) %>%
    select(tf, cell)
  
  # defining cell type combinations
  cell_type_combos <- data.frame(
    atac = atac_ct,
    cell = rna_ct
  )
  
  # go through these cell type combos and ensure cell type naming conventions match those in the ATAC data
  expr_tf_table_mod <- expr_tf_table_new %>%
    left_join(cell_type_combos, by = "cell", relationship = "many-to-many") %>%
    select(tf, atac) %>%
    distinct()
  
  colnames(expr_tf_table_mod)[colnames(expr_tf_table_mod) == "atac"] <- "cell"
  
  # now identifying which cell type-specific TFs identified previously are also being transcribed in the cell type
  expr_prio_tf <- dplyr::intersect(prio_tf_table, expr_tf_table_mod)
  
  return(expr_prio_tf)
}

# read user desired TF expression analysis
tf_expr_req <- args[["tf_expr_analysis"]]

if (tf_expr_req == "atac") {
  # load the necessary libraries
  suppressPackageStartupMessages(library(ape))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(tibble))
  
  # obtain list of TFs to test
  TFs <- tf_table_filt$tf
  TFs <- unique(unlist(strsplit(TFs, "::")))
  
  # defining list of heterodimers
  TF_heterodimers <- unique(tf_table_filt$tf)
  TF_heterodimers <- TF_heterodimers[grepl("::", TF_heterodimers, fixed = T)]
  
  # obtain promoter region definitions and directories necessary for analysis
  prom_thr_up = if (nzchar(args[["prom_th_up"]])) {
    as.integer(args[["prom_th_up"]])
  } else {
    defaults$prom_th_up
  } 
  prom_thr_down = if (nzchar(args[["prom_th_down"]])) {
    as.integer(args[["prom_th_down"]])
  } else {
    defaults$prom_th_down
  }
  gene_annot_path = args[["gencode_dir"]]
  peak_file_dir = paste0(output_dir, "cell_peak.tsv")
  
  tf_expr_results <- confirm_tf_promoter_peaks(TFs, TF_heterodimers, prom_thr_up, prom_thr_down, 
                                               peak_file_dir, gene_annot_path, tf_table_filt)
  
  # accounting for if there are no results found
  if (is.data.frame(tf_expr_results)) {
    # convert this dataframe into binary matrix
    tf_expr_results$value = "atac"
    
    tf_expr_mtx <- pivot_wider(tf_expr_results,
                               names_from = cell,
                               values_from = value,
                               values_fill = "none")
    
    tf_expr_mtx <- column_to_rownames(tf_expr_mtx, var = "tf")
    tf_expr_mtx <- as.matrix(tf_expr_mtx)
    
    # save this matrix
    write.table(tf_expr_mtx, file = paste0(output_dir, "all_TF_expr_results.txt"), row.names = T, quote = F,
                sep = '\t')
  } else {
    message("Note: No TF expression suggested through scATAC-seq analysis. One possible reason for this is a mismatch between genome builds of scATAC-seq and GENCODE data.")
  }
} else if (tf_expr_req == "rna") {
  # load the necessary libraries
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(methods))
  
  # obtain user instructions
  scrna_instruct_dir <- args[["scrna_instruct_dir"]]
  
  # checking to see if scRNA-seq instructions provided
  if (!file.exists(scrna_instruct_dir)) {
    stop("scRNA-seq analysis requested but scRNA-seq instructions spreadsheet not found, please ensure path is correct/provided. Or if scRNA-seq analysis is not desired please set the tf_expr_analysis parameter to 'atac' or do not use this parameter.")
  } 
  
  user_instruct <- read.table(scrna_instruct_dir, header = TRUE)
  
  # getting list of atac cell types requested and rmp cell types, seeing if rmps were not found in some atac cell types, removing them
  rmp_cell_types <- unique(tf_table_filt$cell)
  peak_cell_types <- unique(unlist(str_split(user_instruct$atac_cell_types, ',')))
  missing_cell_types <- setdiff(peak_cell_types, rmp_cell_types)
  
  if (length(missing_cell_types) > 0) { # removing missing cell types from user instructions
    message("Note: TF expression analysis will not be conducted on cell types in which no SNPs co-localizing with risk-mediating peaks were identified: ", paste(missing_cell_types, collapse = ", "), ". Analysis not being conducted may also occur if a cell type name in the atac_cell_types column does not match with the corresponding cell type name provided in the scATAC-seq object (i.e., B cell vs b_cell).")
  }
  
  # obtaining some user input
  user_instruct$one_cell_type_seurat <- as.logical(user_instruct$one_cell_type_seurat)
  quantile_th <- if (nzchar(args[["tf_rna_quantile_th"]])) {
    as.numeric(args[["tf_rna_quantile_th"]])
  } else {
    message("Using default quantile_th value: ", defaults$quantile_th)
    defaults$quantile_th
  }
  
  # creating list of results for scRNA-seq
  scrna_results_list <- list()
  
  # accounting for if scRNA-seq datasets only have one cell type
  if (any(user_instruct$one_cell_type_seurat)) {
    one_ct_entries <- user_instruct[user_instruct$one_cell_type_seurat,]
    
    # create for loop to iterate through each of these entries
    num_one_ct <- nrow(one_ct_entries)
    
    for (i in 1:num_one_ct) {
      current_row <- one_ct_entries[i,]
      scrna_dir <- current_row$scrna_dir
      
      # read in scRNA-seq dataset and associated info
      seurat_name <- load(scrna_dir)
      scrna_dataset <- get(seurat_name)
      atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
      matrix_loc <- current_row$matrix_loc
      
      # filter tf by cell type table and get list of TFs
      tf_table_ct <- tf_table_filt[tf_table_filt$cell %in% atac_cell_types,]
      TFs <- tf_table_ct$tf
      TFs <- unique(unlist(strsplit(TFs, "::")))
      
      # defining list of heterodimers
      TF_heterodimers <- unique(tf_table_ct$tf)
      TF_heterodimers <- TF_heterodimers[grepl("::", TF_heterodimers, fixed = T)]
      
      # conduct TF expression analysis
      percent_df <- gene_percent_expr_ct(scrna_dataset, matrix_loc)
      quant_value <- as.numeric(quant_per_celltype(percent_df, quantile_th))
      expression_df <- tf_expression_analysis_ct(scrna_dataset, TFs, quant_value, atac_cell_types, matrix_loc)
      if (is.data.frame(expression_df)) {
        if (length(TF_heterodimers) > 0) { 
          expr_df_hd <- heterodimers_expr(TF_heterodimers, expression_df)
          tfs_prio_expr <- finalize_tfs_prio_expr_ct(tf_table_ct, expr_df_hd, atac_cell_types)
        } else { # accounting for case where there are no heterodimers
          expression_df <- column_to_rownames(expression_df, var = 'TFs_prio')
          tfs_prio_expr <- finalize_tfs_prio_expr_ct(tf_table_ct, expression_df, atac_cell_types)
        }
        scrna_results_list[[length(scrna_results_list) + 1]] <- tfs_prio_expr
      } else {
        message("Note: No scRNA-seq TF analysis results found for row with scRNA-seq directory '", scrna_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. scRNA-seq dataset may not contain any JASPAR TFs prioritized in the effect-SNP stage or no TFs passed the expression threshold.")
      }
    }
  } 
  
  # now looking at scRNA-seq datasets with more than one cell type
  if (any(!user_instruct$one_cell_type_seurat)) {
    # filter user_instruct for just those entries where one_cell_type_seurat is equal to FALSE
    other_entries <- user_instruct[!user_instruct$one_cell_type_seurat,]
    
    # create for loop to iterate through user instructions and conduct scRNA-seq TF expression analysis 
    num_other <- nrow(other_entries)
    
    for (i in 1:num_other) {
      current_row <- other_entries[i,]
      scrna_dir <- current_row$scrna_dir
      
      # read in scRNA-seq dataset and associated info
      seurat_name <- load(scrna_dir)
      scrna_dataset <- get(seurat_name)
      atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
      rna_cell_types <- unlist(strsplit(current_row$rna_cell_types, split = ","))
      matrix_loc <- current_row$matrix_loc
      
      # filter the scRNA-seq dataset for just those cell types requested by user
      scrna_dataset <- subset(scrna_dataset, idents = rna_cell_types)
      
      # filter tf by cell type table and get list of TFs
      tf_table_ct <- tf_table_filt[tf_table_filt$cell %in% atac_cell_types,]
      TFs <- tf_table_ct$tf
      TFs <- unique(unlist(strsplit(TFs, "::")))
      
      # defining list of heterodimers
      TF_heterodimers <- unique(tf_table_ct$tf)
      TF_heterodimers <- TF_heterodimers[grepl("::", TF_heterodimers, fixed = T)]
      
      # conduct TF expression analysis
      percent_df <- gene_percent_expr(scrna_dataset, matrix_loc)
      quant_vec <- quant_per_celltype(percent_df, quantile_th)
      expression_df <- tf_expression_analysis(scrna_dataset, TFs, quant_vec, matrix_loc)
      if (is.data.frame(expression_df)) {
        if (length(TF_heterodimers) > 0) { 
          expr_df_hd <- heterodimers_expr(TF_heterodimers, expression_df)
          tfs_prio_expr <- finalize_tfs_prio_expr(tf_table_ct, expr_df_hd, atac_cell_types, rna_cell_types)
        } else { # accounting for case where there are no heterodimers
          expression_df <- column_to_rownames(expression_df, var = 'TFs_prio')
          tfs_prio_expr <- finalize_tfs_prio_expr(tf_table_ct, expression_df, atac_cell_types, rna_cell_types)
        }
        scrna_results_list[[length(scrna_results_list) + 1]] <- tfs_prio_expr
      } else {
        message("Note: No scRNA-seq TF analysis results found for row with scRNA-seq directory '", scrna_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. scRNA-seq dataset may not contain any JASPAR TFs prioritized in the effect-SNP stage or no TFs passed the expression threshold.")
      }
    }
  }
  
  # check to see if any results were obtained
  if (length(scrna_results_list) == 0) {
    message('scRNA-seq analysis complete, no results were found over all scRNA-seq datasets.')
  } else {
    # creating one unified dataframe
    all_results <- bind_rows(scrna_results_list) %>% distinct()
    
    # create a binary matrix
    all_results_mtx <- all_results %>%
      mutate(present = 'rna') %>%
      pivot_wider(
        names_from = cell,
        values_from = present,
        values_fill = 'none'
      ) %>%
      column_to_rownames("tf") %>%
      as.matrix()
    
    # save this matrix
    write.table(all_results_mtx, file = paste0(output_dir, "all_TF_expr_results.txt"), row.names = T, quote = F,
                sep = '\t')
  }
} else if (tf_expr_req == "both") {
  # load the necessary libraries
  suppressPackageStartupMessages(library(Seurat))
  suppressPackageStartupMessages(library(tibble))
  suppressPackageStartupMessages(library(methods))
  suppressPackageStartupMessages(library(ape))
  suppressPackageStartupMessages(library(GenomicRanges))
  suppressPackageStartupMessages(library(tibble))
  
  # first conducting RNA-seq analysis
  # obtain user instructions
  scrna_instruct_dir <- args[["scrna_instruct_dir"]]
  
  # checking to see if scRNA-seq instructions provided
  if (!file.exists(scrna_instruct_dir)) {
    stop("scRNA-seq analysis requested but scRNA-seq instructions spreadsheet not found, please ensure path is correct/provided. Or if scRNA-seq analysis is not desired please set the tf_expr_analysis parameter to 'atac' or do not use this parameter.")
  } 
  
  user_instruct <- read.table(scrna_instruct_dir, header = TRUE)
  
  # getting list of atac cell types requested and rmp cell types, seeing if rmps were not found in some atac cell types, removing them
  rmp_cell_types <- unique(tf_table_filt$cell)
  peak_cell_types <- unique(unlist(str_split(user_instruct$atac_cell_types, ',')))
  missing_cell_types <- setdiff(peak_cell_types, rmp_cell_types)
  
  if (length(missing_cell_types) > 0) { # removing missing cell types from user instructions
    message("Note: TF expression analysis will not be conducted on cell types in which no SNPs co-localizing with risk-mediating peaks were identified: ", paste(missing_cell_types, collapse = ", "), ". Analysis not being conducted may also occur if a cell type name in the atac_cell_types column does not match with the corresponding cell type name provided in the scATAC-seq object (i.e., B cell vs b_cell).")
  }
  
  user_instruct$one_cell_type_seurat <- as.logical(user_instruct$one_cell_type_seurat)
  quantile_th <- if (nzchar(args[["tf_rna_quantile_th"]])) {
    as.numeric(args[["tf_rna_quantile_th"]])
  } else {
    message("Using default quantile_th value: ", defaults$quantile_th)
    defaults$quantile_th
  }
  
  # creating list of all TF expression analysis results
  tf_expr_results_list <- list()
  
  # accounting for if scRNA-seq datasets only have one cell type
  if (any(user_instruct$one_cell_type_seurat)) {
    one_ct_entries <- user_instruct[user_instruct$one_cell_type_seurat,]
    
    # create for loop to iterate through each of these entries
    num_one_ct <- nrow(one_ct_entries)
    
    for (i in 1:num_one_ct) {
      current_row <- one_ct_entries[i,]
      scrna_dir <- current_row$scrna_dir
      
      # read in scRNA-seq dataset and associated info
      seurat_name <- load(scrna_dir)
      scrna_dataset <- get(seurat_name)
      atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
      matrix_loc <- current_row$matrix_loc
      
      # filter tf by cell type table and get list of TFs
      tf_table_ct <- tf_table_filt[tf_table_filt$cell %in% atac_cell_types,]
      TFs <- tf_table_ct$tf
      TFs <- unique(unlist(strsplit(TFs, "::")))
      
      # defining list of heterodimers
      TF_heterodimers <- unique(tf_table_ct$tf)
      TF_heterodimers <- TF_heterodimers[grepl("::", TF_heterodimers, fixed = T)]
      
      # conduct TF expression analysis
      percent_df <- gene_percent_expr_ct(scrna_dataset, matrix_loc)
      quant_value <- as.numeric(quant_per_celltype(percent_df, quantile_th))
      expression_df <- tf_expression_analysis_ct(scrna_dataset, TFs, quant_value, atac_cell_types, matrix_loc)
      if (is.data.frame(expression_df)) {
        if (length(TF_heterodimers) > 0) { 
          expr_df_hd <- heterodimers_expr(TF_heterodimers, expression_df)
          tfs_prio_expr <- finalize_tfs_prio_expr_ct(tf_table_ct, expr_df_hd, atac_cell_types)
        } else { # accounting for case where there are no heterodimers
          expression_df <- column_to_rownames(expression_df, var = 'TFs_prio')
          tfs_prio_expr <- finalize_tfs_prio_expr_ct(tf_table_ct, expression_df, atac_cell_types)
        }
        tf_expr_results_list[[length(tf_expr_results_list) + 1]] <- tfs_prio_expr
      } else {
        message("Note: No scRNA-seq TF analysis results found for row with scRNA-seq directory '", scrna_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. scRNA-seq dataset may not contain any JASPAR TFs prioritized in the effect-SNP stage or no TFs passed the expression threshold.")
      }
    }
  } 
  
  # now looking at scRNA-seq datasets with more than one cell type
  if (any(!user_instruct$one_cell_type_seurat)) {
    # filter user_instruct for just those entries where one_cell_type_seurat is equal to FALSE
    other_entries <- user_instruct[!user_instruct$one_cell_type_seurat,]
    
    # create for loop to iterate through user instructions and conduct scRNA-seq TF expression analysis 
    num_other <- nrow(other_entries)
    
    for (i in 1:num_other) {
      current_row <- other_entries[i,]
      scrna_dir <- current_row$scrna_dir
      
      # read in scRNA-seq dataset and associated info
      seurat_name <- load(scrna_dir)
      scrna_dataset <- get(seurat_name)
      atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
      rna_cell_types <- unlist(strsplit(current_row$rna_cell_types, split = ","))
      matrix_loc <- current_row$matrix_loc
      
      # filter the scRNA-seq dataset for just those cell types requested by user
      scrna_dataset <- subset(scrna_dataset, idents = rna_cell_types)
      
      # filter tf by cell type table and get list of TFs
      tf_table_ct <- tf_table_filt[tf_table_filt$cell %in% atac_cell_types,]
      TFs <- tf_table_ct$tf
      TFs <- unique(unlist(strsplit(TFs, "::")))
      
      # defining list of heterodimers
      TF_heterodimers <- unique(tf_table_ct$tf)
      TF_heterodimers <- TF_heterodimers[grepl("::", TF_heterodimers, fixed = T)]
      
      # conduct TF expression analysis
      percent_df <- gene_percent_expr(scrna_dataset, matrix_loc)
      quant_vec <- quant_per_celltype(percent_df, quantile_th)
      expression_df <- tf_expression_analysis(scrna_dataset, TFs, quant_vec, matrix_loc)
      if (is.data.frame(expression_df)) {
        if (length(TF_heterodimers) > 0) { 
          expr_df_hd <- heterodimers_expr(TF_heterodimers, expression_df)
          tfs_prio_expr <- finalize_tfs_prio_expr(tf_table_ct, expr_df_hd, atac_cell_types, rna_cell_types)
        } else { # accounting for case where there are no heterodimers
          expression_df <- column_to_rownames(expression_df, var = 'TFs_prio')
          tfs_prio_expr <- finalize_tfs_prio_expr(tf_table_ct, expression_df, atac_cell_types, rna_cell_types)
        }
        tf_expr_results_list[[length(tf_expr_results_list) + 1]] <- tfs_prio_expr
      } else {
        message("Note: No scRNA-seq TF analysis results found for row with scRNA-seq directory '", scrna_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. scRNA-seq dataset may not contain any JASPAR TFs prioritized in the effect-SNP stage or no TFs passed the expression threshold.")
      }
    }
  }
  
  # now conducting ATAC-seq analysis 
  # obtain list of TFs to test
  TFs <- tf_table_filt$tf
  TFs <- unique(unlist(strsplit(TFs, "::")))
  
  # defining list of heterodimers
  TF_heterodimers <- unique(tf_table_filt$tf)
  TF_heterodimers <- TF_heterodimers[grepl("::", TF_heterodimers, fixed = T)]
  
  # obtain promoter region definitions and directories necessary for analysis
  prom_thr_up = if (nzchar(args[["prom_th_up"]])) {
    as.integer(args[["prom_th_up"]])
  } else {
    defaults$prom_th_up
  } 
  prom_thr_down = if (nzchar(args[["prom_th_down"]])) {
    as.integer(args[["prom_th_down"]])
  } else {
    defaults$prom_th_down
  }
  gene_annot_path = args[["gencode_dir"]]
  peak_file_dir = paste0(output_dir, "cell_peak.tsv")
  
  # predict TF expression using ATAC-seq
  tf_expr_results <- confirm_tf_promoter_peaks(TFs, TF_heterodimers, prom_thr_up, prom_thr_down, 
                                               peak_file_dir, gene_annot_path, tf_table_filt)
  
  # bringing both analyses together
  if (is.data.frame(tf_expr_results) & length(tf_expr_results_list) > 0) { # both rna-seq and atac-seq results are available
    tf_expr_results$conf <- 'atac'
    
    # creating one unified dataframe
    all_results <- bind_rows(tf_expr_results_list) %>% distinct()
    all_results$conf <- "rna"
    all_results <- bind_rows(all_results, tf_expr_results)
    
    # creating matrix
    all_results_val <- all_results %>%
      group_by(tf, cell) %>%
      summarise(
        support_type = paste(sort(unique(conf)), collapse = ","),
        .groups = "drop"
      ) %>%
      mutate(
        value = case_when(
          support_type == "atac" ~ "atac",
          support_type == "rna" ~ "rna",
          support_type == "atac,rna" ~ "both"
        )
      ) %>%
      select(tf, cell, value) %>% 
      distinct()
    
    all_results_mtx <- all_results_val %>%
      pivot_wider(
        names_from = cell,
        values_from = value,
        values_fill = "none"
      ) %>%
      column_to_rownames(var = 'tf') %>%
      as.matrix()
    
    # save this matrix
    write.table(all_results_mtx, file = paste0(output_dir, "all_TF_expr_results.txt"), row.names = T, quote = F,
                sep = '\t')
    
    message('TF expression analysis complete!')
  } else if (is.data.frame(tf_expr_results)) { # just atac-seq results available
    # convert this dataframe into matrix
    tf_expr_results$value = "atac"
    
    tf_expr_mtx <- pivot_wider(tf_expr_results,
                               names_from = cell,
                               values_from = value,
                               values_fill = "none")
    
    tf_expr_mtx <- column_to_rownames(tf_expr_mtx, var = "tf")
    tf_expr_mtx <- as.matrix(tf_expr_mtx)
    
    # save this matrix
    write.table(tf_expr_mtx, file = paste0(output_dir, "all_TF_expr_results.txt"), row.names = T, quote = F,
                sep = '\t')
    
    print('TF expression analysis complete! Results found for scATAC-seq analysis. No results found for scRNA-seq analysis.')
  } else if (length(tf_expr_results_list) > 0) { # just rna-seq results available
    # creating one unified dataframe
    all_results <- bind_rows(tf_expr_results_list) %>% distinct()
    
    # create a matrix
    all_results_mtx <- all_results %>%
      mutate(present = 'rna') %>%
      pivot_wider(
        names_from = cell,
        values_from = present,
        values_fill = 'none'
      ) %>%
      column_to_rownames("tf") %>%
      as.matrix()
    
    # save this matrix
    write.table(all_results_mtx, file = paste0(output_dir, "all_TF_expr_results.txt"), row.names = T, quote = F,
                sep = '\t')
    
    message('TF expression analysis complete! No results found for scATAC-seq analysis. Results found for scRNA-seq analysis.')
  } else { # no results found
    message('TF expression analysis complete, no results were found over both the scRNA-seq and scATAC-seq analyses.')
  }
} 

message("TF expression analysis complete!")