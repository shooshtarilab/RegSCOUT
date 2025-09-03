suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(xlsx))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(Repitools))
suppressPackageStartupMessages(library(reshape2))
suppressPackageStartupMessages(library(dplyr))

#Setting command line arguments
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Defining default parameter values
defaults <- list(
  peak_th = 0.1, 
  cell_count_th = 3
)

#Loading the scATAC-seq object 
pbmc_dir = args[["seurat_obj"]]
pbmc_name = load(pbmc_dir)
pbmc = get(pbmc_name)
rm(pbmc_name)

#Removing the cell types with less than a certain number of cells in them
cell_count_th <- if (nzchar(args[["cell_count_th"]])) {
  as.numeric(args[["cell_count_th"]])
} else {
  message("Using default cell_count_th value: ", defaults$cell_count_th)
  defaults$cell_count_th
}

for (cell in levels(x = pbmc)){
  if(ncol(subset(x = pbmc, idents = cell)) < cell_count_th){
    print(paste0(cell, " had less than ", cell_count_th, " cells"))
    pbmc= subset(x = pbmc, idents = cell, invert = TRUE)
  }else{
    print(paste0(cell, " had equal or more than ", cell_count_th, " cells"))
  }
}
print(paste0("Cell types with less than ", cell_count_th, " cells removed!"))

#Getting the list of cell types
cell_types = unique(Idents(pbmc))

#Initializing cell by peak matrix
loci_cluster_matrix = matrix(0, ncol = length(cell_types), nrow = nrow(pbmc))
row.names(loci_cluster_matrix) = row.names(pbmc)
colnames(loci_cluster_matrix) = cell_types

#Extracting counts matrix from scATAC-seq object
pbmc_counts = pbmc@assays$peaks$counts

#Extracting cell type specific peaks
print("Extracting cell-type specific regions:")
peak_th = if (nzchar(args[["peak_th"]])) {
  as.numeric(args[["peak_th"]])
} else {
  message("Using default peak_th value: ", defaults$peak_th)
  defaults$peak_th
}

for (cell in cell_types){
  cell_index = Idents(pbmc) == cell
  cell_counts = pbmc_counts[,cell_index]
  
  if (inherits(pbmc_counts, "dgCMatrix")) { # account for sparse matrices
    cell_counts@x = rep(1,length(cell_counts@x))
  } else {
    cell_counts[which(cell_counts>0)] = 1
  }
  cell_counts[which(cell_counts>0)] = 1
  cell_counts_perc = rowSums(cell_counts)/ncol(cell_counts)
  cell_peaks = names(cell_counts_perc)[cell_counts_perc > peak_th]
  
  loci_cluster_matrix[cell_peaks,cell] = 1
}

#Filtering the cell by peak matrix to only include the peaks
#found in at least one cell type
loci_cluster_matrix = loci_cluster_matrix[rowSums(loci_cluster_matrix)>0,]

#Turning the cell by peak matrix to a dataframe showing the cell type
#of each peak
loci_cluster_df <- melt(loci_cluster_matrix)
loci_cluster_df = loci_cluster_df[loci_cluster_df$value > 0,]
loci_cluster_df = loci_cluster_df[,1:2]
colnames(loci_cluster_df) = c("Peak","Cell_Type")

loci_cluster_df <- loci_cluster_df %>%
  group_by(Peak) %>%
  summarize(Cell_Type = paste(unique(Cell_Type), collapse = ","))

peak_sep_frame = str_split(loci_cluster_df$Peak, pattern = "-",simplify = TRUE)

final_peak_cell_df = as.data.frame(matrix(0, nrow = nrow(loci_cluster_df),
                                          ncol = 4))
colnames(final_peak_cell_df) = c("chr","start","end","cell sub-types")

final_peak_cell_df$chr = peak_sep_frame[,1]
final_peak_cell_df$start = peak_sep_frame[,2]
final_peak_cell_df$end = peak_sep_frame[,3]
final_peak_cell_df$`cell sub-types` = loci_cluster_df$Cell_Type

#Getting the output directory
output_file_main = args[["output_dir"]]

#Saving the final table in the output directory
output_file = paste0(output_file_main, "cell_peak.xlsx")
write.xlsx(final_peak_cell_df, file = output_file, col.names = TRUE,
           row.names = FALSE)
print('Finished extracting cell type-specific open chromatin regions!')
