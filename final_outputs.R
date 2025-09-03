library(tibble)
library(stringr)
library(tidyr)
library(dplyr)
library(readxl)
library(GenomicRanges)
library(ape)
library(writexl)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# read output directory
output_dir <- args[["output_dir"]]

# defining default parameter values
defaults <- list(
  gene_ppa_th = 0.05
)

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
if (!is.null(user_finemap)) {
  fine_loci_head <- read.table(ci_gwas_dir, sep = '\t', header = T)
  
  # reading in loci info file and adding its information to fine_loci_head
  loci_info <- read.table(paste0(output_dir, 'loci_info.txt'), sep = '\t', header = T)
  
  # bring locus information into fine_loci_head
  fine_loci_head <- merge(fine_loci_head, loci_info, by = "chunk")
} else {
  fine_loci_head <- read.table(ci_gwas_dir, sep = '\t', header = T)
}


## preparing tf data for effect_snps located on rmps (associating SNPs with TFs)
risk_regions_ratio = read.table(paste0(output_dir, 'risk_regions_ratio.txt'), header = TRUE)
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

# adding locus regions
for (i in 1:nrow(fine_loci_head)) {
  fine_loci_head$locus_region[i] <- paste0(fine_loci_head$locus_chr[i], "-", fine_loci_head$locus_start[i], "-", fine_loci_head$locus_end[i])
}

# adding in locus region and lead SNP ID
for (i in 1:nrow(final_output)) {
  current_locus <- final_output$locus[i]
  match_index <- which(fine_loci_head$chunk == current_locus)
  match_index <- match_index[1]
  
  if (length(match_index) > 0) {
    final_output$locus_region[i] <- fine_loci_head$locus_region[match_index]
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
risk_regions_ppa = read.table(paste0(output_dir, "risk_regions_ppa.txt"), header = TRUE)

# adding in rmp information
for (i in 1:nrow(final_output)) {
  current_rmp <- final_output$rmp[i]
  match_index <- which(risk_regions_ppa$region == current_rmp)
  if (length(match_index) > 0) {
    final_output$rmp_ppa[i] <- risk_regions_ppa$sumPPA[match_index]
    final_output$cell_type[i] <- risk_regions_ppa$cell[match_index]
  }
}

# adding promoter-linked rmps genes information
rmp_promoter = read.table(paste0(output_dir, "direct_rmp_gene_overlaps.txt"),header=TRUE)
rmp_promoter <- rmp_promoter[,c("RMP", "Gene")] %>% distinct()

rmp_gene <- rmp_promoter %>% group_by(RMP) %>% summarise(gene = paste0(Gene, collapse = ", "))

for (i in 1:nrow(final_output)) {
  current_rmp <- final_output$rmp[i]
  match_index <- which(rmp_gene$RMP == current_rmp)
  if (length(match_index) > 0) {
    final_output$directly_mapped_gene[i] <- rmp_gene$gene[match_index]
  }
}

# one row per cell type
final_output <- final_output %>% separate_rows(cell_type, sep = ",")

# reading in cicero information
cicero = read.table(paste0(output_dir, "cic_peak_interact_gene.txt"), header = TRUE)

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

# bringing in Hi-C results, if Hi-C analysis was requested
hic_eqtl_req <- args[["hic_eqtl_analysis"]]

if (!is.null(hic_eqtl_req)) {
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
}

# adding in histone mark information if histone mark analysis was requested by user
hist_mark_req <- args[["histone_mark_analysis"]]

if (tolower(hist_mark_req) == "y") {
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
if (!is.null(hic_eqtl_req)) {
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
}

write.table(final_output, file = paste0(output_dir, "final_table.txt") ,row.names = FALSE, quote = FALSE, sep = "\t")

# now creating the prioritized table
final_table <- final_output

# creating new columns
final_table$gene_score <- NA
final_table$gene <- NA

# finding all genes identified in locus
for (i in 1:nrow(final_table)) {
  final_table$gene[i] <- paste0(final_table$directly_mapped_gene[i], ", ", final_table$cicero_gene[i], ", ", final_table$hic_gene[i], ", ", final_table$eqtl_gene[i])
}

# creating a separate row for each gene that could be prioritized by cicero, direct overlap, hic, or eqtl
final_table <- final_table %>% separate_rows(gene, sep = ", ")
final_table <- final_table[!final_table$gene %in% "NA",] # remove the rows with NA in the gene column
final_table <- final_table %>% distinct()

# calculating gene ppa values
gene_ppa_df <- final_table[,c('cell_type', 'gene', 'rmp', 'rmp_ppa')] %>% distinct()
gene_ppa_df <- gene_ppa_df %>%
  group_by(cell_type, gene) %>%
  summarise(gene_ppa = sum(rmp_ppa, na.rm = T), .groups = 'drop')

final_table <- final_table %>%
  left_join(gene_ppa_df, by = c("cell_type", "gene"))

# gene prioritization score
scoring_system <- c(
  "directly_mapped_gene" = 2,
  "cicero_gene" = 1,
  "hic_gene" = 1,
  "eqtl_gene" = 1
)

for (i in seq_len(nrow(final_table))) {
  gene <- final_table$gene[i]
  score <- 0
  
  for (col in names(scoring_system)) {
    col_val <- final_table[[col]][i]
    
    # Split by ", " and trim spaces
    genes_in_col <- trimws(unlist(strsplit(col_val, ",\\s*")))
    
    if (gene %in% genes_in_col) {
      score <- score + scoring_system[[col]]
    }
  }
  
  final_table$gene_score[i] <- score
}

# now establishing TF priority scores
# first making it so that each row only has one TF listed
final_table_new <- final_table %>%
  separate_rows(tf, sep = ", ") %>%
  mutate(tf = str_trim(tf))

# establishing separate columns for TF direction and TF name
final_table_new$tf_dir <- sub("^[^\\(]+", "", final_table_new$tf)
final_table_new$tf <- gsub("\\(.*?\\)", "", final_table_new$tf)

# reading in TF confirmation results if TF confirmation was conducted
TF_conf_req <- args[["tf_expr_analysis"]]

if (TF_conf_req == "both") {
  # reading in TF confirmation analysis results
  tf_confirmed <- read.table(paste0(output_dir, "all_TF_expr_results.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # converting this matrix into a dataframe
  tf_confirmed_df <- tf_confirmed %>%
    rownames_to_column(var = "tf") %>%
    pivot_longer(
      cols = -tf,
      names_to = "cell_type",
      values_to = "evidence"
    ) %>%
    filter(evidence > 0)
  
  # now introducing two columns tf_peak and tf_rna for ease of integration into final table
  tf_confirmed_df <- tf_confirmed_df %>%
    mutate(
      tf_peak = case_when(
        evidence == 1 ~ 1,
        evidence == 2 ~ 0,
        evidence == 3 ~ 1
      ),
      tf_rna = case_when(
        evidence == 1 ~ 0,
        evidence == 2 ~ 1,
        evidence == 3 ~ 1
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
  
} else if (TF_conf_req == 'atac' | TF_conf_req == 'rna') {
  # reading in TF confirmation analysis results
  tf_confirmed <- read.table(paste0(output_dir, "all_TF_expr_results.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # converting this matrix into a dataframe and naming column based on whether atac or rna conf was requested
  if (TF_conf_req == 'atac') {
    tf_confirmed_df <- tf_confirmed %>%
      rownames_to_column(var = "tf") %>%
      pivot_longer(
        cols = -tf,
        names_to = "cell_type",
        values_to = "tf_peak"
      ) %>%
      filter(tf_peak > 0)
  } else if (TF_conf_req == 'rna') {
    tf_confirmed_df <- tf_confirmed %>%
      rownames_to_column(var = "tf") %>%
      pivot_longer(
        cols = -tf,
        names_to = "cell_type",
        values_to = "tf_rna"
      ) %>%
      filter(tf_rna > 0)
  }
  
  # bringing in that TF expression information
  final_table_new <- final_table_new %>%
    left_join(
      tf_confirmed_df,
      by = c("tf", "cell_type")
    )
}

write.table(final_table_new, file = paste0(output_dir, "complete_table.txt") ,row.names = FALSE, quote = FALSE, sep = "\t")


# now performing gene filtering, prioritization
if (!is.null(hic_eqtl_req)) {
  gene_ppa_th <- if (!is.null(args[["gene_ppa_th"]])) as.numeric(args[["gene_ppa_th"]]) else defaults$gene_ppa_th
  prioritized_table <- final_table_new %>%
    filter(gene_score > 1) %>%
    filter(gene_ppa > gene_ppa_th) %>%
    select(-c(hic_prom, cicero_prom)) %>%
    distinct()
} else { # if Hi-C or eQTL analysis was not requested, don't filter by gene score
  prioritized_table <- final_table_new %>%
    filter(gene_ppa > gene_ppa_th) %>%
    distinct()
}

# printing a summary of findings for this dataset
summary_stat <- c("SNPs" = length(unique(prioritized_table$effect_snp)), 
                  "Genes" = length(unique(prioritized_table$gene)), 
                  "Cell Types" = length(unique(prioritized_table$cell_type)), 
                  "TFs" = length(unique(prioritized_table$tf)))
print("Number of unique SNPs, genes, cell types, and TFs in prioritized table:")
print(summary_stat)

# save prioritized table
write.table(prioritized_table, paste0(output_dir, "prioritized_table.txt"),row.names = FALSE, quote = FALSE, sep = "\t")







print('Final tables created!')
