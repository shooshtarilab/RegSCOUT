suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(circlize))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# read output directory
output_dir <- args[["output_dir"]]

# defining default parameter values
defaults <- list(
  gene_sum_ppa_th = 0.05,
  gene_score_th = 2,
  tf_score_th = 0
)

# Remove gencode file generated in Peak_Gene_Integration, used in TF expression analysis and HIC
invisible(suppressWarnings(file.remove(paste0(output_dir, "gene_tss_granges.rds"))))

# starting to create the final dataframe
final_output <- data.frame(matrix(ncol = 17, nrow = 1))

colnames(final_output) <- c("locus", "locus_region",	"lead_snp",	"lead_ppa",	
                            "effect_snp",	"effect_ppa",	"tf","rmp",	"rmp_hm_label", "rmp_ppa", 
                            "cicero_prom", "hic_prom", "directly_mapped_gene",	"cicero_gene",	
                            "hic_gene",	"eqtl_gene", "cell_type")

# read in ci_gwas_dir
ci_gwas_dir = args[["ci_gwas_dir"]]

# read in if fine-mapping was desired
user_finemap <- args[["finemap"]]

## assigning the locus number to each ci_snp, read in CI (credible interval) SNP file
if (user_finemap == "Y") {
  fine_loci_head <- read.table(ci_gwas_dir, sep = '\t', header = T)
  
  # reading in loci info file and adding its information to fine_loci_head
  loci_info <- read.table(paste0(output_dir, 'loci_info.txt'), sep = '\t', header = T)
  
  # bring locus information into fine_loci_head
  fine_loci_head <- merge(fine_loci_head, loci_info, by = "chunk")
} else {
  fine_loci_head <- read.table(ci_gwas_dir, sep = '\t', header = T)
}

## preparing tf data for effect_snps located on rmps (associating SNPs with TFs)
risk_regions_ratio <- read.table(paste0(output_dir, 'risk_regions_ratio.txt'),header=TRUE)

risk_regions_ratio <- risk_regions_ratio %>%
  mutate(
    tf = str_extract(TFSNP, "^.+(?=-[^-]+$)"), # accounts for TFs e.g., NKX6-2 with dash in name
    snp = str_extract(TFSNP, "[^\\-]+$") 
  )

## defining impact of SNP on tf based on log_lik_ratio
risk_regions_ratio$direction <- NA

for (i in 1:nrow(risk_regions_ratio)) {
  if (risk_regions_ratio$log_lik_ratio[i] > 0) {
    risk_regions_ratio$direction[i] <- "(+)"
  } else {
    risk_regions_ratio$direction[i] <- "(-)"
  }
}

risk_regions_ratio$tf_direction <- paste0(risk_regions_ratio$tf, risk_regions_ratio$direction)

# adding snp and tf information to the final_output data frame
final_output <- final_output[rep(1, nrow(risk_regions_ratio)),]
final_output$effect_snp <- risk_regions_ratio$snp
final_output$tf <- risk_regions_ratio$tf_direction
final_output <- final_output %>% distinct()

# one row per effect SNP, group all TFs affected tgt
final_output <- final_output %>%
  group_by(effect_snp) %>%
  mutate(tf = paste(unique(tf), collapse = ", ")) %>%
  distinct()

# adding in effect SNP PPA and locus ids
for (i in 1:nrow(final_output)) {
  current_snp <- final_output$effect_snp[i]
  match_index <- which(fine_loci_head$id == current_snp)
  if (length(match_index) > 0) {
    final_output$effect_ppa[i] <- fine_loci_head$PPA[match_index]
    final_output$locus[i] <- fine_loci_head$chunk[match_index]
  }
}

# adding locus regions if that information is available
if (all(c("locus_chr", "locus_start", "locus_end") %in% colnames(fine_loci_head))) {
  for (i in 1:nrow(fine_loci_head)) {
    fine_loci_head$locus_region[i] <- paste0(fine_loci_head$locus_chr[i], "-", fine_loci_head$locus_start[i], "-", fine_loci_head$locus_end[i])
  }
  
  for (i in 1:nrow(final_output)) {
    current_locus <- final_output$locus[i]
    match_index <- which(fine_loci_head$chunk == current_locus)
    match_index <- match_index[1]
    
    if (length(match_index) > 0) {
      final_output$locus_region[i] <- fine_loci_head$locus_region[match_index]
    }
  }
} else {
  message('No information on locus regions found in fine-mapped data, this information will not be included in final tables. It is recommended to provide this information if available.')
  final_output$locus_region <- NULL
}

# adding in lead snp and lead snp ppa if that information is available
if ('lead_snp' %in% colnames(fine_loci_head)) {
  for (i in 1:nrow(final_output)) {
    current_locus <- final_output$locus[i]
    match_index <- which(fine_loci_head$chunk == current_locus)
    match_index <- match_index[1]
    
    if (length(match_index) > 0) {
      final_output$lead_snp[i] <- fine_loci_head$lead_snp[match_index]
    }
  }
  
  # adding in lead SNP PPA
  for (i in 1:nrow(final_output)) {
    current_lead <- final_output$lead_snp[i]
    match_index <- which(fine_loci_head$id == current_lead)
    if (length(match_index) > 0) {
      final_output$lead_ppa[i] <- fine_loci_head$PPA[match_index]
    }
  }
} else {
  message("No information on lead SNPs found in fine-mapped data, this information will not be included in final tables. It is recommended to provide this information if available.")
  final_output$lead_snp <- NULL
  final_output$lead_ppa <- NULL
}

final_output <- final_output[order(final_output$locus),]

# adding rmp information
snp_rmp <- risk_regions_ratio[,c('snp','region')] %>% distinct()
snp_rmp <- snp_rmp %>% group_by(snp) %>% summarise(rmp = paste0(unique(region), collapse = ", "))

for (i in 1:nrow(final_output)) {
  current_snp <- final_output$effect_snp[i]
  match_index <- which(snp_rmp$snp == current_snp)
  if (length(match_index) > 0) {
    final_output$rmp[i] <- snp_rmp$rmp[match_index]
  }
}

# for the rare case where one SNP overlaps with two RMPs (it may fall at the edge of both)
final_output <- final_output %>% separate_rows(rmp, sep = ", ")

# loading in rmps and their PPAs (associated with the E-SNPs that overlap them)
risk_regions_ppa <- read.table(paste0(output_dir, "risk_regions_ppa.txt"), header=TRUE)

# adding in rmp information
for (i in 1:nrow(final_output)) {
  current_rmp <- final_output$rmp[i]
  match_index <- which(risk_regions_ppa$region == current_rmp)
  if (length(match_index) > 0) {
    final_output$rmp_ppa[i] <- risk_regions_ppa$sumPPA[match_index]
    final_output$cell_type[i] <- risk_regions_ppa$cell[match_index]
  }
}

# adding promoter-linked rmps genes information, if any were found
direct_ovlp_res_dir <- paste0(output_dir, "direct_rmp_gene_overlaps.txt")

if (file.exists(direct_ovlp_res_dir)) {
  rmp_promoter <- read.table(direct_ovlp_res_dir, header = TRUE, sep = '\t')
  
  rmp_promoter <- rmp_promoter[,c("RMP", "Gene")] %>% distinct()
  
  rmp_gene <- rmp_promoter %>% group_by(RMP) %>% summarise(gene = paste0(Gene, collapse = ", "))
  
  for (i in 1:nrow(final_output)) {
    current_rmp <- final_output$rmp[i]
    match_index <- which(rmp_gene$RMP == current_rmp)
    if (length(match_index) > 0) {
      final_output$directly_mapped_gene[i] <- rmp_gene$gene[match_index]
    }
  }
} else {
  final_output$directly_mapped_gene <- NULL
}

# one row per cell type
final_output <- final_output %>% separate_rows(cell_type, sep = ",")

# reading in cicero information, if results were found
cicero_dir <- paste0(output_dir, "cic_peak_interact_gene.txt")

if (file.exists(cicero_dir)) {
  cicero <- read.table(cicero_dir, header = TRUE, sep = '\t')
  
  # renaming certain columns
  colnames(cicero)[colnames(cicero) == "Peak1"] <- 'rmp'
  colnames(cicero)[colnames(cicero) == "Peak2"] <- 'prom_peak'
  colnames(cicero)[colnames(cicero) == "Cell_Type"] <- 'cell_type'
  
  # filter this dataframe to add cicero_prom and cicero_gene information
  cicero_gene_df <- cicero[,c('rmp', 'prom_peak', 'cell_type', 'Gene')] %>% distinct()
  
  # group genes that had same promoter peak
  cicero_gene_df <- cicero_gene_df %>%
    group_by(rmp, prom_peak, cell_type) %>%
    summarise(
      Gene = unique(Gene) %>%
        sort() %>%
        str_c(collapse = ", "),
      .groups = "drop"
    )
  
  # adding the cicero information
  final_output <- final_output %>%
    left_join(
      cicero_gene_df,
      by = c("rmp", "cell_type"), 
      relationship = "many-to-many"
    )
  
  final_output$cicero_prom <- final_output$prom_peak 
  final_output$cicero_gene <- final_output$Gene 
  final_output$prom_peak <- NULL
  final_output$Gene <- NULL
} else {
  final_output$cicero_prom <- NULL
  final_output$cicero_gene <- NULL
}

# bringing in Hi-C results, if Hi-C analysis was requested
hic_results_dir <- paste0(output_dir, "all_hic_results.txt")

if (file.exists(hic_results_dir)) {
  hic_results <- read.table(paste0(output_dir, "all_hic_results.txt"), sep = '\t', header = T)
  
  # removing unnecessary column for this step 
  hic_results$oeRegion <- NULL
  hic_results <- hic_results %>% distinct()
  
  # renaming some columns
  colnames(hic_results)[colnames(hic_results) == "rmpRegion"] <- "rmp"
  colnames(hic_results)[colnames(hic_results) == "cell"] <- "cell_type"
  
  # group genes that had same promoter peak
  hic_gene_df <- hic_results %>%
    group_by(rmp, promRegion, cell_type) %>%
    summarise(
      gene = unique(gene) %>%
        sort() %>%
        str_c(collapse = ", "),
      .groups = "drop"
    )
  
  # adding the hic information
  final_output <- final_output %>%
    left_join(
      hic_gene_df,
      by = c("rmp", "cell_type"), 
      relationship = "many-to-many"
    )
  
  final_output$hic_prom <- final_output$promRegion 
  final_output$hic_gene <- final_output$gene
  final_output$promRegion <- NULL
  final_output$gene <- NULL
} else {
  final_output$hic_prom <- NULL
  final_output$hic_gene <- NULL
}

# adding in histone mark information if histone mark analysis was requested by user
his_mark_results_dir <- paste0(output_dir, "all_histone_mark_results.txt")

if (file.exists(his_mark_results_dir)) {
  hist_mark_results <- read.table(paste0(output_dir, 'all_histone_mark_results.txt'), sep = '\t', header = T)
  
  # remove unnecessary columns for this step 
  hist_mark_results$version <- NULL 
  
  # renaming certain columns
  colnames(hist_mark_results)[colnames(hist_mark_results) == 'cell'] <- 'cell_type'
  
  # adding the histone mark information
  final_output <- final_output %>%
    left_join(
      hist_mark_results,
      by = c("rmp", "cell_type")
    )
  
  final_output$rmp_hm_label <- final_output$state
  final_output$state <- NULL
} else {
  final_output$rmp_hm_label <- NULL # remove column if histone mark analysis was not requested
}

# bringing in eQTL results 
eqtl_results_dir <- paste0(output_dir, 'all_eqtl_results.txt')

if (file.exists(eqtl_results_dir)) {
  eqtl_results <- read.table(paste0(output_dir, 'all_eqtl_results.txt'), sep = '\t', header = T)
  
  # remove and rename some columns for this step
  eqtl_for_table <- eqtl_results
  eqtl_for_table$rmp <- NULL
  eqtl_for_table <- eqtl_for_table %>% distinct()
  colnames(eqtl_for_table)[colnames(eqtl_for_table) == "snp"] <- "effect_snp"
  colnames(eqtl_for_table)[colnames(eqtl_for_table) == "cell"] <- "cell_type"
  
  # group genes tgt if they are prioritized with the same SNP in the same cell type
  eqtl_for_table <- eqtl_for_table %>%
    group_by(effect_snp, cell_type) %>%
    summarise(
      gene = unique(gene) %>%
        sort() %>%
        str_c(collapse = ", "),
      .groups = "drop"
    )
  
  # adding the eqtl information
  final_output <- final_output %>%
    left_join(
      eqtl_for_table,
      by = c("effect_snp", "cell_type")
    )
  
  final_output$eqtl_gene <- final_output$gene
  final_output$gene <- NULL
} else {
  final_output$eqtl_gene <- NULL
}

# now creating the prioritized table
final_table <- final_output

# if no lines of evidence for genes were found, do not run rest of code
gene_columns <- c("directly_mapped_gene",	"cicero_gene", "hic_gene", "eqtl_gene")
cols_present <- intersect(gene_columns, colnames(final_table))

if (length(cols_present) == 0) {
  stop('Final table created, no lines of evidence for any risk-mediating genes were found.')
}

# creating new columns
final_table$gene_score <- NA
final_table$gene <- NA

# finding all genes identified in locus
for (i in 1:nrow(final_table)) {
  # add values in those columns to gene_vec
  gene_vec <- as.vector(final_table[i,cols_present,drop = F])
  gene_char <- paste0(gene_vec, collapse = ", ")
  
  # add this to gene column of final table
  final_table$gene[i] <- gene_char
}

# creating a separate row for each gene that could be prioritized by cicero, direct overlap, hic, or eqtl
final_table <- final_table %>% separate_rows(gene, sep = ", ")
final_table <- final_table[!final_table$gene %in% "NA",] # remove the rows with NA in the gene column
final_table <- final_table %>% distinct()

# calculating gene ppa values
gene_sum_ppa_df <- final_table[,c('cell_type', 'gene', 'rmp', 'rmp_ppa')] %>% distinct()
gene_sum_ppa_df <- suppressMessages(gene_sum_ppa_df %>%
  group_by(cell_type, gene) %>%
  summarise(gene_sum_ppa = sum(rmp_ppa, na.rm = T), .groups = 'drop'))

final_table <- final_table %>%
  left_join(gene_sum_ppa_df, by = c("cell_type", "gene"))

# gene prioritization score, filtering scoring system by what columns are present in final table
scoring_system <- c(
  "directly_mapped_gene" = 2,
  "cicero_gene" = 1,
  "hic_gene" = 1,
  "eqtl_gene" = 1
)

scoring_system <- scoring_system[cols_present]

# now establishing gene score
for (i in seq_len(nrow(final_table))) {
  gene <- final_table$gene[i]
  score <- 0
  
  for (col in names(scoring_system)) {
    col_val <- final_table[[col]][i]
    
    # Split by ", " and trim whitespace
    genes_in_col <- trimws(unlist(strsplit(col_val, ",\\s*")))
    
    if (gene %in% genes_in_col) {
      score <- score + scoring_system[[col]]
    }
  }
  
  final_table$gene_score[i] <- score
}

# now establishing TF priority scores
tf_expr_results_dir <- paste0(output_dir, "all_TF_expr_results.txt")

final_table_new <- final_table

if (file.exists(tf_expr_results_dir)) {
  # first making it so that each row only has one TF listed
  final_table_new <- final_table_new %>%
    separate_rows(tf, sep = ", ") %>%
    mutate(tf = str_trim(tf))
  
  # establishing separate columns for TF direction and TF name
  final_table_new$tf_dir <- sub("^[^\\(]+", "", final_table_new$tf)
  final_table_new$tf <- gsub("\\(.*?\\)", "", final_table_new$tf)
  
  # read in tf confirmation analysis results
  tf_confirmed <- read.table(paste0(output_dir, "all_TF_expr_results.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # get all unique values possible in the matrix
  type_tf_conf <- unique(unlist(tf_confirmed))
  
  # reading in TF confirmation results to create TF heatmaps and incorporate results in final_table_new
  if (all(c("rna", "atac") %in% type_tf_conf) | "both" %in% type_tf_conf) { # if both atac and rna confirmation were done
    ## creating a TF by cell type heatmap with tf expression information
    # creating a TF by cell type matrix from final_table_new
    tf_cell_df <- risk_regions_ratio[,c('cell', 'tf')] %>% distinct()
    
    tf_cell_matrix <- tf_cell_df %>%
      mutate(value = 1) %>%                
      pivot_wider(
        names_from = cell, 
        values_from = value,
        values_fill = list(value = 0)
      ) %>% 
      column_to_rownames(var = 'tf') %>%
      as.matrix()
    
    # adding missing cell types into tf_confirmed_prep if any missing and matching column orders
    cell_types <- colnames(tf_cell_matrix)
    missing_cts <- setdiff(cell_types, colnames(tf_confirmed))
    tf_confirmed_prep <- tf_confirmed
    
    if (length(missing_cts) > 0) {
      for (ct in missing_cts) {
        tf_confirmed_prep[[ct]] <- "none"
      }
    }
    
    # reorder columns to match cell_types order
    tf_confirmed_prep <- tf_confirmed_prep[, cell_types]
    
    # adding in missing TFs into tf_confirmed_prep if any missing
    missing_tfs <- setdiff(rownames(tf_cell_matrix), rownames(tf_confirmed_prep))
    
    if (length(missing_tfs) > 0) {
      tfs_to_add <- as.data.frame(
        matrix("none", 
               nrow = length(missing_tfs), 
               ncol = ncol(tf_confirmed_prep), 
               dimnames = list(missing_tfs, colnames(tf_confirmed_prep)))
      )
      
      tf_confirmed_prep <- rbind(tf_confirmed_prep, tfs_to_add)
    }
    
    # ensure rows are ordered alphabetically
    tf_cell_matrix <- tf_cell_matrix[order(rownames(tf_cell_matrix)),]
    tf_confirmed_prep <- tf_confirmed_prep[order(rownames(tf_confirmed_prep)),]
    tf_confirmed_prep <- as.matrix(tf_confirmed_prep)
    
    # turning character values into numerical values in tf_confirmed_prep
    tf_confirmed_prep[tf_confirmed_prep == 'none'] <- 0
    tf_confirmed_prep[tf_confirmed_prep == 'atac'] <- 1
    tf_confirmed_prep[tf_confirmed_prep == 'rna'] <- 2
    tf_confirmed_prep[tf_confirmed_prep == 'both'] <- 3
    rownames_mtx <- rownames(tf_confirmed_prep)
    tf_confirmed_prep <- apply(tf_confirmed_prep, 2, as.numeric)
    rownames(tf_confirmed_prep) <- rownames_mtx
    
    # creating one matrix that reflects all information in the two matrices created
    tf_heatmap_mtx <- tf_confirmed_prep + tf_cell_matrix
    
    # create heatmap
    f1 = colorRamp2(seq(0, 4, length = 5), c("#EEEEEE", "#ffffcc", "#a1dab4", "#41b6c4", "#225ea8"))
    legend_levels <- c('no evidence', 'snp impact', 'snp impact + atac', 'snp impact + rna', 'snp impact + atac + rna')
    
    heatmap_TFs <- Heatmap(
      tf_heatmap_mtx, name = "TF heatmap", col = f1, 
      column_title = "Cell Type", row_title = "Risk-Mediating TF",
      column_title_gp = grid::gpar(fontsize = 14),
      row_title_gp = grid::gpar(fontsize = 14),
      row_names_gp = grid::gpar(fontsize = 11),
      column_names_gp = grid::gpar(fontsize = 13),
      rect_gp = gpar(col= "#84878a"),
      column_title_side = "bottom",
      heatmap_legend_param = list(
        title = "Evidence",
        at = c(0, 1, 2, 3, 4),
        labels = legend_levels,
        color_bar = "discrete",
        legend_width = unit(16, "cm")
      ),
      width = unit(1 * ncol(tf_heatmap_mtx), "cm"),
      height = unit(1 * nrow(tf_heatmap_mtx), "cm")
    ) 
    
    output_tf_file = paste0(output_dir,"cell_tf.svg")
    svg(output_tf_file, width = (ncol(tf_heatmap_mtx) + 10)/2.54, height = (nrow(tf_heatmap_mtx) + 10)/2.54)
    print(heatmap_TFs)
    invisible(dev.off())
    
    ## bringing TF expression results into final_table_new
    # converting the TF conf. matrix into a dataframe
    tf_confirmed_df <- tf_confirmed %>%
      rownames_to_column(var = "tf") %>%
      pivot_longer(
        cols = -tf,
        names_to = "cell_type",
        values_to = "evidence"
      ) %>%
      filter(evidence != 'none')
    
    # now introducing two columns tf_peak and tf_rna for ease of integration into final table
    tf_confirmed_df <- tf_confirmed_df %>%
      mutate(
        tf_peak = case_when(
          evidence == 'atac' ~ 1,
          evidence == 'rna' ~ 0,
          evidence == 'both' ~ 1
        ),
        tf_rna = case_when(
          evidence == 'atac' ~ 0,
          evidence == 'rna' ~ 1,
          evidence == 'both' ~ 1
        )
      )
    
    tf_confirmed_df$evidence <- NULL
    
    # bringing in that TF expression information
    final_table_new <- final_table_new %>%
      left_join(
        tf_confirmed_df,
        by = c("tf", "cell_type")
      )
    
    # calculating a TF score
    final_table_new$tf_score <- rowSums(final_table_new[, c("tf_peak", "tf_rna")], na.rm = TRUE)
  } else {
    ## creating a TF by cell type heatmap with tf expression information
    # creating a TF by cell type matrix from final_table_new
    tf_cell_df <- risk_regions_ratio[,c('cell', 'tf')] %>% distinct()
    
    tf_cell_matrix <- tf_cell_df %>%
      mutate(value = 1) %>%                
      pivot_wider(
        names_from = cell, 
        values_from = value,
        values_fill = list(value = 0)
      ) %>% 
      column_to_rownames(var = 'tf') %>%
      as.matrix()
    
    # adding missing cell types into tf_confirmed_prep if any missing and matching column orders
    cell_types <- colnames(tf_cell_matrix)
    missing_cts <- setdiff(cell_types, colnames(tf_confirmed))
    tf_confirmed_prep <- tf_confirmed
    
    if (length(missing_cts) > 0) {
      for (ct in missing_cts) {
        tf_confirmed_prep[[ct]] <- "none"
      }
    }
    
    # reorder columns to match cell_types order
    tf_confirmed_prep <- tf_confirmed_prep[, cell_types]
    
    # adding in missing TFs into tf_confirmed_prep if any missing
    missing_tfs <- setdiff(rownames(tf_cell_matrix), rownames(tf_confirmed_prep))
    
    if (length(missing_tfs) > 0) {
      tfs_to_add <- as.data.frame(
        matrix("none", 
               nrow = length(missing_tfs), 
               ncol = ncol(tf_confirmed_prep), 
               dimnames = list(missing_tfs, colnames(tf_confirmed_prep)))
      )
      
      tf_confirmed_prep <- rbind(tf_confirmed_prep, tfs_to_add)
    }
    
    # ensure rows are ordered alphabetically
    tf_cell_matrix <- tf_cell_matrix[order(rownames(tf_cell_matrix)),]
    tf_confirmed_prep <- tf_confirmed_prep[order(rownames(tf_confirmed_prep)),]
    tf_confirmed_prep <- as.matrix(tf_confirmed_prep)
    
    # turning character values into numerical values in tf_confirmed_prep
    tf_confirmed_prep[tf_confirmed_prep == 'none'] <- 0
    tf_confirmed_prep[tf_confirmed_prep == 'atac' | tf_confirmed_prep == 'rna'] <- 1
    rownames_mtx <- rownames(tf_confirmed_prep)
    tf_confirmed_prep <- apply(tf_confirmed_prep, 2, as.numeric)
    rownames(tf_confirmed_prep) <- rownames_mtx
    
    # creating one matrix that reflects all information in the two matrices created
    tf_heatmap_mtx <- tf_confirmed_prep + tf_cell_matrix
    
    # create heatmap
    if ('atac' %in% type_tf_conf) {
      f1 = colorRamp2(seq(0, 2, length = 3), c("#EEEEEE", "#ffffcc", "#a1dab4"))
      legend_levels <- c('no evidence', 'snp impact', 'snp impact + atac')
    } else {
      f1 = colorRamp2(seq(0, 2, length = 3), c("#EEEEEE", "#ffffcc", "#41b6c4"))
      legend_levels <- c('no evidence', 'snp impact', 'snp impact + rna')
    }
    
    heatmap_TFs <- Heatmap(
      tf_heatmap_mtx, name = "TF heatmap", col = f1, 
      column_title = "Cell Type", row_title = "Risk-Mediating TF",
      column_title_gp = grid::gpar(fontsize = 14),
      row_title_gp = grid::gpar(fontsize = 14),
      row_names_gp = grid::gpar(fontsize = 11),
      column_names_gp = grid::gpar(fontsize = 13),
      rect_gp = gpar(col= "#84878a"),
      column_title_side = "bottom",
      heatmap_legend_param = list(
        title = "Evidence",
        at = c(0, 1, 2),
        labels = legend_levels,
        color_bar = "discrete",
        legend_width = unit(16, "cm")
      ),
      width = unit(1 * ncol(tf_heatmap_mtx), "cm"),
      height = unit(1 * nrow(tf_heatmap_mtx), "cm")
    ) 
    
    output_tf_file = paste0(output_dir,"cell_tf.svg")
    svg(output_tf_file, width = (ncol(tf_heatmap_mtx) + 10)/2.54, height = (nrow(tf_heatmap_mtx) + 10)/2.54)
    print(heatmap_TFs)
    invisible(dev.off())
    
    ## bringing TF expression results into final_table_new
    # converting this matrix into a dataframe and naming column based on whether atac or rna conf was requested
    if ('atac' %in% type_tf_conf) {
      tf_confirmed_df <- tf_confirmed %>%
        rownames_to_column(var = "tf") %>%
        pivot_longer(
          cols = -tf,
          names_to = "cell_type",
          values_to = "tf_peak"
        ) %>%
        filter(tf_peak != 'none')
      tf_confirmed_df$tf_peak[tf_confirmed_df$tf_peak == 'atac'] <- 1
    } else if ('rna' %in% type_tf_conf) {
      tf_confirmed_df <- tf_confirmed %>%
        rownames_to_column(var = "tf") %>%
        pivot_longer(
          cols = -tf,
          names_to = "cell_type",
          values_to = "tf_rna"
        ) %>%
        filter(tf_rna != 'none')
      tf_confirmed_df$tf_rna[tf_confirmed_df$tf_rna == 'rna'] <- 1
    }
    
    # bringing in that TF expression information
    final_table_new <- final_table_new %>%
      left_join(
        tf_confirmed_df,
        by = c("tf", "cell_type")
      )
  }
}

write.table(final_table_new, file = paste0(output_dir, "complete_table.txt"), row.names = FALSE, quote = FALSE, sep = "\t")

# now performing gene filtering, prioritization
# first determining if filtering by tf score is desired
tf_score_th  <- if (is.null(args[["tf_score_th"]])) {
  message("Using default tf_score_th value: ", defaults$tf_score_th)
  defaults$tf_score_th
} else if (!nzchar(args[["tf_score_th"]])) {
  message("Using default tf_score_th value: ", defaults$tf_score_th)
  defaults$tf_score_th
} else {
  as.integer(args[["tf_score_th"]])
}

gene_sum_ppa_th  <- if (is.null(args[["gene_sum_ppa_th"]])) {
  message("Using default gene_sum_ppa_th value: ", defaults$gene_sum_ppa_th)
  defaults$gene_sum_ppa_th
} else if (!nzchar(args[["gene_sum_ppa_th"]])) {
  message("Using default gene_sum_ppa_th value: ", defaults$gene_sum_ppa_th)
  defaults$gene_sum_ppa_th
} else {
  as.numeric(args[["gene_sum_ppa_th"]])
}

gene_score_th  <- if (is.null(args[["gene_score_th"]])) {
  message("Using default gene_score_th value: ", defaults$gene_score_th)
  defaults$gene_score_th
} else if (!nzchar(args[["gene_score_th"]])) {
  message("Using default gene_score_th value: ", defaults$gene_score_th)
  defaults$gene_score_th
} else {
  as.numeric(args[["gene_score_th"]])
}

if (tf_score_th > 0 & file.exists(tf_expr_results_dir)) {
  if (all(c('rna', 'atac') %in% type_tf_conf)) {
    prioritized_table <- final_table_new %>%
      filter(gene_score >= gene_score_th) %>%
      filter(tf_score >= tf_score_th) %>%
      filter(gene_sum_ppa >= gene_sum_ppa_th) %>%
      distinct()
  } else if ('atac' %in% type_tf_conf) {
    prioritized_table <- final_table_new %>%
      filter(gene_score >= gene_score_th) %>%
      filter(tf_peak >= tf_score_th) %>%
      filter(gene_sum_ppa >= gene_sum_ppa_th) %>%
      distinct()
  } else if ('rna' %in% type_tf_conf) {
    prioritized_table <- final_table_new %>%
      filter(gene_score >= gene_score_th) %>%
      filter(tf_rna >= tf_score_th) %>%
      filter(gene_sum_ppa >= gene_sum_ppa_th) %>%
      distinct()
  }
} else {
  prioritized_table <- final_table_new %>%
    filter(gene_score >= gene_score_th) %>%
    filter(gene_sum_ppa >= gene_sum_ppa_th) %>%
    distinct()
}

if (nrow(prioritized_table) == 0) {
  stop("Final tables created, no prioritized genes found with gene_score >= ", gene_score_th, ", gene_sum_ppa_th >= ", gene_sum_ppa_th, ", and tf_score_th >= ", tf_score_th, ".")
}

if (file.exists(hic_results_dir) & file.exists(cicero_dir)) { 
  prioritized_table <- prioritized_table %>%
    select(-c(hic_prom, cicero_prom)) %>%
    distinct()
} else if (file.exists(hic_results_dir)){ 
  prioritized_table <- prioritized_table %>%
    select(-c(hic_prom)) %>%
    distinct()
} else if (file.exists(cicero_dir)) {
  prioritized_table <- prioritized_table %>%
    select(-c(cicero_prom)) %>%
    distinct()
}

# creating a prioritized gene by cell type heat map, starting with creating a data frame for the plot
gene_ct_df <- prioritized_table[,c('cell_type','gene','gene_score')] %>% distinct()
gene_ct_df <- suppressMessages(gene_ct_df %>% # obtaining max gene_score value for each cell type - gene pair
  group_by(cell_type, gene) %>%
  summarise(gene_score = max(gene_score)) %>% 
  ungroup())

# turning gene_ct_df into a matrix
gene_cell_mtx <- gene_ct_df %>%
  pivot_wider(
    names_from = cell_type, 
    values_from = gene_score,
    values_fill = list(gene_score = 0)
  ) %>% 
  column_to_rownames(var = 'gene') %>%
  as.matrix()

# getting all possible values of gene_score
break_points <- sort(unique(as.vector(gene_cell_mtx)))

# creating a colour palette with length equal to number of breakpoints and creating circlize function
color_palette <- colorRampPalette(c("#EEEEEE", "#d4b9da", "#c994c7", "#df65b0", "#dd1c77", "#980043"))(length(break_points))
f2 <- colorRamp2(break_points, color_palette)

# creating heatmap
legend_levels <- as.character(break_points)

heatmap_genes <- Heatmap(
  gene_cell_mtx, name = "Gene heatmap", col = f2, 
  column_title = "Cell Type", row_title = "Prioritized Gene",
  column_title_gp = grid::gpar(fontsize = 14),
  row_title_gp = grid::gpar(fontsize = 14),
  row_names_gp = grid::gpar(fontsize = 11),
  column_names_gp = grid::gpar(fontsize = 13),
  rect_gp = gpar(col= "#84878a"),
  column_title_side = "bottom",
  heatmap_legend_param = list(
    title = "Confidence Score",
    at = break_points,
    labels = legend_levels,
    color_bar = "discrete",
    legend_width = unit(16, "cm")
  ),
  width = unit(1 * ncol(gene_cell_mtx), "cm"),
  height = unit(1 * nrow(gene_cell_mtx), "cm")
) 

output_gene_file = paste0(output_dir,"cell_gene.svg")
svg(output_gene_file, width = (ncol(gene_cell_mtx) + 10)/2.54, height = (nrow(gene_cell_mtx) + 10)/2.54)
print(heatmap_genes)
invisible(dev.off())

# printing a summary of findings for this dataset
summary_stat <- c("Effect-SNPs" = length(unique(prioritized_table$effect_snp)), 
                  "Genes" = length(unique(prioritized_table$gene)), 
                  "Cell Types" = length(unique(prioritized_table$cell_type)), 
                  "TFs" = length(unique(unlist(strsplit(prioritized_table$tf, ", ")))))
message("Number of unique effect-SNPs, genes, cell types, and TFs in prioritized table:")
message(paste(summary_stat, collapse = "   "))

# save prioritized table
write.table(prioritized_table, file = paste0(output_dir, "prioritized_table.txt"), quote = FALSE, sep = "\t")

# creating a condensed version of the prioritized table
if (file.exists(tf_expr_results_dir)) {
  summary_table <- prioritized_table %>%
    mutate(tf_dir = gsub("[()]", "", tf_dir))
} else {
  # first making it so that each row only has one TF listed
  prioritized_table <- prioritized_table %>%
    separate_rows(tf, sep = ", ") %>%
    mutate(tf = str_trim(tf))
  
  # establishing separate columns for TF direction and TF name
  prioritized_table$tf_dir <- sub("^[^\\(]+", "", prioritized_table$tf)
  prioritized_table$tf <- gsub("\\(.*?\\)", "", prioritized_table$tf)
  
  summary_table <- prioritized_table %>%
    mutate(tf_dir = gsub("[()]", "", tf_dir))
}

# establish a new column with tfs and a summary of their scores/directions
if (file.exists(tf_expr_results_dir)) {
  if (all(c("rna", "atac") %in% type_tf_conf) | "both" %in% type_tf_conf) {
    summary_table <- summary_table %>%
      mutate(`tfs (direction, score)` = paste0(tf, " (", tf_dir, ", ", tf_score, ")"))
  } else if ("rna" %in% type_tf_conf) {
    summary_table$tf_rna[is.na(summary_table$tf_rna)] <- 0 # such that there are no NA values in the new column
    summary_table <- summary_table %>%
      mutate(`tfs (direction, score)` = paste0(tf, " (", tf_dir, ", ", tf_rna, ")"))
  } else if ("atac" %in% type_tf_conf) {
    summary_table$tf_peak[is.na(summary_table$tf_peak)] <- 0 # such that there are no NA values in the new column
    summary_table <- summary_table %>%
      mutate(`tfs (direction, score)` = paste0(tf, " (", tf_dir, ", ", tf_peak, ")"))
  }
} else {
  summary_table <- summary_table %>%
    mutate(`tfs (direction)` = paste0(tf, " (", tf_dir, ")"))
}

# remove certain columns
columns_to_remove <- c('tf_score', 'tf_dir', 'tf', 'directly_mapped_gene', 'cicero_gene', 'hic_gene',
            'eqtl_gene', 'tf_peak', 'tf_rna')

summary_table <- summary_table %>%
  select(-any_of(columns_to_remove)) %>% 
  distinct()

# summarizing some stats relevant to genes before putting them into summary table
summary_table$gene_sum_ppa <- round(summary_table$gene_sum_ppa, 2)
summary_table <- summary_table %>%
  group_by(effect_snp, rmp, cell_type, gene) %>%
  mutate(gene_score = max(gene_score)) %>%
  ungroup() %>%
  distinct()

# establish a new column with genes and a summary of their scores/ppas
summary_table <- summary_table %>%
  mutate(`genes (score, sum_ppa)` = paste0(gene, " (", gene_score, ", ", gene_sum_ppa, ")"))

# remove certain columns
summary_table <- summary_table %>%
  select(-c(gene_score, gene_sum_ppa, gene)) %>% 
  distinct()

# grouping TFs and genes based on effect SNP, RMP, and cell type
if (file.exists(tf_expr_results_dir)) {
  table_results_combined <- summary_table %>%
    group_by(effect_snp, rmp, cell_type) %>%
    mutate(
      `tfs (direction, score)` = paste(unique(`tfs (direction, score)`), collapse = ", "),
      `genes (score, sum_ppa)` = paste(unique(`genes (score, sum_ppa)`), collapse = ", ")
    ) %>%
    ungroup() %>%
    distinct()
} else {
  table_results_combined <- summary_table %>%
    group_by(effect_snp, rmp, cell_type) %>%
    mutate(
      `tfs (direction)` = paste(unique(`tfs (direction)`), collapse = ", "),
      `genes (score, sum_ppa)` = paste(unique(`genes (score, sum_ppa)`), collapse = ", ")
    ) %>%
    ungroup() %>%
    distinct()
}

# cleaning up some numerical values
table_results_combined <- table_results_combined %>%
  mutate(effect_ppa = round(effect_ppa, 4),
         lead_ppa = round(lead_ppa, 4))

# write this file out
write.table(table_results_combined, file = paste0(output_dir, "prioritized_table_condensed.txt"), quote = FALSE, sep = "\t")