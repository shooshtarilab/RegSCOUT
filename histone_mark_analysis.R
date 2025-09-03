suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(GenomicRange))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(purrr))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# read output directory
output_dir = args[["output_dir"]]

# read in user instructions
hist_mark_instruct_dir <- args[["hist_mark_instruct_dir"]]

# loading in user instructions
user_instruct <- read.table(hist_mark_instruct_dir, header = TRUE)

# organizing user instructions
atac_cell_types <- user_instruct$atac_cell_types
hist_mark_cell_types <- user_instruct$chromHMM_cell_types
hist_dirs <- user_instruct$chromHMM_dir

# initiating a list for chromHMM data
hist_data_list <- setNames(vector("list", length(hist_mark_cell_types)), hist_mark_cell_types)

col_names <- c("chr", "start", "end", "state") # filtering for required columns

for (i in seq_along(hist_dirs)) {
  hist_data_list[[i]] <- read.table(hist_dirs[[i]], header = T, sep = "\t")
  hist_data_list[[i]] <- hist_data_list[[i]][,col_names]
}

# creating GRanges objects for chromHMM data
hist_granges <- list()

for (cell in names(hist_data_list)) {
    hist_granges[[cell]] <- GRanges(seqnames = hist_data_list[[cell]]$chr, ranges = IRanges(start = hist_data_list[[cell]]$start, end = hist_data_list[[cell]]$end))
    hist_granges[[cell]]$state = hist_data_list[[cell]]$state
}

# reading in risk-mediating peak regions
peak_cluster_matrix = read.delim(paste0(output_dir, "peak_cluster_matrix.txt"), header = TRUE, sep = "\t")

# preparing cell subtype specific RMP data
result_matrices <- list()

rmp_cell_types <- unique(unlist(strsplit(atac_cell_types, split = ","))) # may only have chromHMM data for a subset of atac cell types

for (cell_type in rmp_cell_types) {
  sub_matrix <- peak_cluster_matrix[, cell_type, drop = FALSE]
  filtered_matrix <- sub_matrix[sub_matrix[[1]] == 1, , drop = FALSE]
  result_matrices[[cell_type]] <- filtered_matrix
}

for (i in seq_along(result_matrices)) {
  result_matrices[[i]]$rmp <- rownames(result_matrices[[i]])
  rownames(result_matrices[[i]]) <- NULL
} 

# creating GRanges objects for RMP data
rmp_granges <- list()

for (cell_type in names(result_matrices)) {
  rmp_granges[[cell_type]] <- GRanges(seqnames = sapply(result_matrices[[cell_type]]$rmp, function(x) str_split(x, "-")[[1]][1]), 
                                      ranges = IRanges(start = as.numeric(sapply(result_matrices[[cell_type]]$rmp, function(x) str_split(x, "-")[[1]][2])), 
                                                       end = as.numeric(sapply(result_matrices[[cell_type]]$rmp, function(x) str_split(x, "-")[[1]][3]))))
}

# making cell type pairings
ct_pairs <- data.frame(
  atac = atac_cell_types,
  hist = hist_mark_cell_types
)

ct_pairs <- ct_pairs %>%
  separate_rows(atac, sep = ",", convert = FALSE)

match_cell_types <- ct_pairs %>%
  pull(hist) %>%
  setNames(ct_pairs$atac)

# adding ChromHMM states to RMP data
rmp_labels_list <- list()

for (cell_type in names(match_cell_types)) {
  rmp_grange <- rmp_granges[[cell_type]]
  hist_grange <- hist_granges[[match_cell_types[[cell_type]]]]
  
  rmp_hist_index <- findOverlaps(rmp_grange, hist_grange)
  
  query_hits <- queryHits(rmp_hist_index)
  subject_hits <- subjectHits(rmp_hist_index)
  
  query_granges <- rmp_grange[query_hits]
  subject_granges <- hist_grange[subject_hits]
  state_values <- mcols(subject_granges)$state
  
  rmp_labels <- data.frame(
    rmp = as.character(query_granges),
    hist = as.character(subject_granges),
    state = as.character(state_values)
  )
  
  rmp_labels$rmp <- gsub(":", "-", rmp_labels$rmp)
  rmp_labels$hist <- NULL
  
  # make one row per RMP, group all unique states tgt
  rmp_labels <- rmp_labels %>%
    group_by(rmp) %>%
    summarise(
      state = paste(unique(state), collapse = ", "),
      .groups = "drop"
    )
  
  rmp_hist <- rmp_labels
  
  # determining the chromHMM model used
  rmp_hist$version <- paste0(length(levels(as.factor(hist_granges[[match_cell_types[[cell_type]]]]$state))), '-state') 
  
  rmp_labels_list[[cell_type]] <- rmp_hist
}

# adding cell type labels to each dataframe
rmp_labels_list <- imap(rmp_labels_list, ~ mutate(.x, cell = .y))

# bind all results together
all_results <- bind_rows(rmp_labels_list)

# save this dataframe
write.table(all_results, file = paste0(output_dir, "all_histone_mark_results.txt"), row.names = F, quote = F,
            sep = '\t')


print('Histone mark analysis complete!')
