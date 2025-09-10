suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(ComplexHeatmap))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(circlize))
suppressPackageStartupMessages(library(rtracklayer))

#Setting command line arguments
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Defining default parameter values
defaults <- list(
  prom_th_up = 2000,
  prom_th_down = 2000
)

#Getting the working directory
output_file_main = args[["output_dir"]]

#Getting the cell by peak table
cell_peak_file = paste0(output_file_main,"cell_peak.tsv")
cell_peak = read.delim(cell_peak_file, header = TRUE)
peaks_cell_types = unique(unlist(str_split(cell_peak$cell_types, ",")))

#Getting the file of effect SNPs and loading them
eff_snp_file = paste0(output_file_main,"Ci_effect_SNPs.txt")
eff_snp = read.table(file = eff_snp_file, header = TRUE, sep = '\t')
eff_snp_filt = eff_snp[,c('SNP', 'CHR', 'Pos', 'ppa')]
eff_snp_filt = eff_snp_filt %>% distinct()

#Finding the sum of PPA of peaks
peak_list = paste(cell_peak$chr,cell_peak$start,cell_peak$end,
                  sep = "-")
peak_granges = StringToGRanges(peak_list)
eff_snp_list = paste(eff_snp_filt$CHR,eff_snp_filt$Pos,eff_snp_filt$Pos,
                     sep = "-")
eff_snp_granges = StringToGRanges(eff_snp_list)

peak_eff_overlaps = findOverlaps(peak_granges,eff_snp_granges)

if (length(peak_eff_overlaps) == 0) {
  stop("No co-localization between effect-SNPs and cell type-specific open chromatin regions found. One possible reason for this is a mismatch between genome builds of scATAC-seq and GWAS datasets.")
}

peak_ppa_frame = as.data.frame(matrix(0,nrow = length(peak_eff_overlaps),
                                      ncol = 3))
colnames(peak_ppa_frame) = c("region","PPA","cell")

peak_ppa_frame$region = peak_list[queryHits(peak_eff_overlaps)]
peak_ppa_frame$cell = cell_peak$cell_types[queryHits(peak_eff_overlaps)]
peak_ppa_frame$PPA = eff_snp_filt$ppa[subjectHits(peak_eff_overlaps)]

peak_ppa_frame = peak_ppa_frame %>%
  group_by(region, cell) %>%
  mutate(sumPPA = sum(PPA)) %>%
  ungroup()

peak_ppa_frame = peak_ppa_frame %>% select(-PPA)
peak_ppa_frame = unique(peak_ppa_frame)

rmp_cell_types = unique(unlist(str_split(peak_ppa_frame$cell, ",")))
cell_type_diff = setdiff(peaks_cell_types, rmp_cell_types)

# check to see whether there were any cases where an RMP was not identified in a cell type
if (length(cell_type_diff) > 0) {
  message("Note: No risk-mediating peaks were identified in these cell types: ", paste(cell_type_diff, collapse = ", "), ". One possible reason for this is a mismatch between genome builds of scATAC-seq and GWAS datasets.")
}

peak_ppa_frame_dir = paste0(output_file_main,"risk_regions_ppa.txt")
write.table(as.data.frame(peak_ppa_frame), file = peak_ppa_frame_dir, sep = "\t", quote = FALSE, row.names = FALSE)

peak_cluster_matrix = peak_ppa_frame[,c("region","cell")]
peak_cluster_matrix <- peak_cluster_matrix %>%
  separate_rows(cell, sep = ",")
peak_cluster_matrix <- peak_cluster_matrix %>%
  mutate(value = 1) %>%
  spread(key = cell, value = value, fill = 0)
peak_cluster_matrix = tibble::column_to_rownames(peak_cluster_matrix, var = "region")
peak_cluster_matrix = as.matrix(peak_cluster_matrix)

peak_cluster_matrix_file = paste0(output_file_main,"peak_cluster_matrix.txt")
write.table(peak_cluster_matrix, peak_cluster_matrix_file, sep = "\t", row.names = TRUE, quote = FALSE)

f1 = colorRamp2(seq(0, 1, length = 2), c("#EEEEEE", "red"))


heatmap_peaks <- Heatmap(
  peak_cluster_matrix, name = "Percentage accessibility", col = f1, 
  column_title = "Cell Type", row_title = "Risk-Mediating Peak",
  column_title_gp = grid::gpar(fontsize = 14),
  row_title_gp = grid::gpar(fontsize = 14),
  row_names_gp = grid::gpar(fontsize = 11),
  column_names_gp = grid::gpar(fontsize = 13),
  rect_gp = gpar(col= "#84878a"),
  column_title_side = "bottom",
  heatmap_legend_param = list(title = "Presence", at = c(0, 1), 
                              labels = c("0", "1"), 
                              color_bar = "vertical", 
                              legend_width = unit(16, "cm")),
  width = unit(1 * ncol(peak_cluster_matrix), "cm"),
  height = unit(1 * nrow(peak_cluster_matrix), "cm")
  )


# guideline bar based on sumPPA
risk_regions = peak_ppa_frame[,c("region", "sumPPA")]
rownames(risk_regions) = NULL
risk_regions = tibble::column_to_rownames(risk_regions, var = "region")
risk_regions = as.matrix(risk_regions)

# sorting the rownames of the risk_regions matrix to match the rownames of the peak_cluster_matrix
roworder <- order(match(rownames(risk_regions), rownames(peak_cluster_matrix)))
risk_regions <- risk_regions[roworder, , drop = FALSE]

f2 = colorRamp2(seq(0, max(risk_regions[,1], na.rm = TRUE), length = 2), c("#EEEEEE","#54278f"))

# create a separate heatmap for the ppa
heatmap_ppa <- Heatmap(
  risk_regions, name = "SumPPA", col = f2, 
  heatmap_legend_param = list(title = "SumPPA", 
                              at = c(0, max(risk_regions, na.rm = TRUE)), 
                              labels = c("0", round(max(risk_regions, na.rm = TRUE))), 
                              color_bar = "vertical", 
                              legend_width = unit(16, "cm")),
  row_names_gp = grid::gpar(fontsize = 11),
  column_names_gp = grid::gpar(fontsize = 13),
  # Add cell_fun argument to display values inside the heatmap cells
  layer_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(sprintf("%.3f", risk_regions[i, j]), x, y, gp = grid::gpar(fontsize = 8))
  },
  width = unit(1 * ncol(risk_regions), 'cm'),
  height = unit(1 * nrow(risk_regions), 'cm')
)

# combine the two heatmaps and save it 
output_peak_file = paste0(output_file_main,"cell_peak.svg")
svg(output_peak_file, width = ((ncol(peak_cluster_matrix) + ncol(risk_regions)) + 10) / 2.54, height = (nrow(peak_cluster_matrix) + 10)/ 2.54)
combined_heatmaps <- HeatmapList(heatmap_peaks + heatmap_ppa)
print(combined_heatmaps)
invisible(dev.off())

#Changing the cell by peak matrix so that each line has only one cell type

#Getting the peaks and Effect SNPs as GRanges and overlapping them to find
#SNP-affected peaks and mapping TFs to those peaks 
cell_peak_grg = StringToGRanges(paste(cell_peak$chr,cell_peak$start,
                                      cell_peak$end,sep = "-"))
eff_snp_grg = StringToGRanges(paste(eff_snp$CHR,eff_snp$Pos,
                                    eff_snp$Pos,sep = "-"))

peak_snp_overlap = findOverlaps(cell_peak_grg, eff_snp_grg)


cell_peak_filt = as.data.frame(matrix(0,nrow = length(peak_snp_overlap),
                                      ncol=(ncol(cell_peak)+3)))

colnames(cell_peak_filt) = c(colnames(cell_peak),"SNP","TF","log_lik_ratio")

cell_peak_filt[,colnames(cell_peak)] = cell_peak[queryHits(peak_snp_overlap),]
cell_peak_filt$SNP = eff_snp$SNP[subjectHits(peak_snp_overlap)]
cell_peak_filt$TF = eff_snp$TF[subjectHits(peak_snp_overlap)]
cell_peak_filt$log_lik_ratio = abs(eff_snp$log_like_ratio[subjectHits(peak_snp_overlap)])

cell_tf_snp = as.data.frame(matrix(0, nrow = nrow(cell_peak_filt),
                                   ncol = 4))
colnames(cell_tf_snp) = c("region","cell","TFSNP","log_like_ratio")
cell_tf_snp$region = paste(cell_peak_filt$chr,
                           cell_peak_filt$start,
                           cell_peak_filt$end,
                           sep = "-")
cell_tf_snp$cell = cell_peak_filt$cell_types
cell_tf_snp$TFSNP = paste0(cell_peak_filt$TF,"-",cell_peak_filt$SNP)
cell_tf_snp$log_like_ratio = eff_snp$log_like_ratio[subjectHits(peak_snp_overlap)]

cell_tf_snp = cell_tf_snp %>% separate_rows(cell, sep = ",")

peak_ratio_frame_dir = paste0(output_file_main,"risk_regions_ratio.txt")
write.table(as.data.frame(cell_tf_snp), file = peak_ratio_frame_dir,
           quote = FALSE, sep="\t", row.names = FALSE)

tf_cluster_matrix = cell_tf_snp[,c("TFSNP","cell")]
tf_cluster_matrix <- tf_cluster_matrix %>% distinct() 
tf_cluster_matrix <- tf_cluster_matrix %>%
  mutate(value = 1) %>%
  spread(key = cell, value = value, fill = 0)
tf_cluster_matrix = tibble::column_to_rownames(tf_cluster_matrix, var = "TFSNP")
tf_cluster_matrix = as.matrix(tf_cluster_matrix)

TF_cluster_matrix_file = paste0(output_file_main,"TF_cluster_matrix.txt")
write.table(tf_cluster_matrix, TF_cluster_matrix_file, sep = "\t", quote = FALSE, row.names = TRUE)


f1 = colorRamp2(seq(0, 1, length = 2), c("#EEEEEE", "blue"))


heatmap_TFs <- Heatmap(
  tf_cluster_matrix, name = "Percentage accessibility", col = f1, 
  column_title = "Cell Type", row_title = "TF-SNP Pair",
  column_title_gp = grid::gpar(fontsize = 14),
  row_title_gp = grid::gpar(fontsize = 14),
  row_names_gp = grid::gpar(fontsize = 11),
  column_names_gp = grid::gpar(fontsize = 13),
  rect_gp = gpar(col= "#84878a"),
  column_title_side = "bottom",
  heatmap_legend_param = list(title = "Presence", at = c(0, 1), 
                              labels = c("0", "1"), 
                              color_bar = "vertical", 
                              legend_width = unit(16, "cm")),
  width = unit(1 * ncol(tf_cluster_matrix), "cm"),
  height = unit(1 * nrow(tf_cluster_matrix), "cm")
)  


# guideline bar based on sum log like ratio
risk_tfs = cell_tf_snp[,c("TFSNP", "log_lik_ratio")]
risk_tfs = unique(risk_tfs)
rownames(risk_tfs) = NULL
risk_tfs <- risk_tfs %>%
  group_by(TFSNP) %>%
  summarize(log_lik_ratio = max(abs(log_lik_ratio)), .groups = 'drop')
risk_tfs = tibble::column_to_rownames(risk_tfs, var = "TFSNP")
risk_tfs = as.matrix(risk_tfs)

# sorting the rownames of the risk_regions matrix to match the rownames of the peak_cluster_matrix
roworder <- order(match(rownames(risk_tfs), rownames(tf_cluster_matrix)))
risk_tfs <- risk_tfs[roworder, , drop = FALSE]

f2 = colorRamp2(seq(0, max(risk_tfs[,1], na.rm = TRUE), length = 2), c("#EEEEEE","#e6550d"))

# create a separate heatmap for the ppa
heatmap_ppa <- Heatmap(
  risk_tfs, name = "Log_like_ratio", col = f2, 
  heatmap_legend_param = list(title = "Log_like_ratio", 
                              at = c(0, max(risk_tfs, na.rm = TRUE)), 
                              labels = c("0", round(max(risk_tfs, na.rm = TRUE))), 
                              color_bar = "vertical", 
                              legend_width = unit(16, "cm")),
  row_names_gp = grid::gpar(fontsize = 11),
  column_names_gp = grid::gpar(fontsize = 13),
  # Add cell_fun argument to display values inside the heatmap cells
  layer_fun = function(j, i, x, y, width, height, fill) {
    grid::grid.text(sprintf("%.3f", risk_tfs[i, j]), x, y, gp = grid::gpar(fontsize = 8))
  },
  width = unit(1 * ncol(risk_tfs), 'cm'),
  height = unit(1 * nrow(risk_tfs), 'cm')
)

# combine the two heatmaps and save it 
output_tf_file = paste0(output_file_main,"cell_snp_tf.svg")
svg(output_tf_file,
    width = (ncol(risk_tfs) + ncol(tf_cluster_matrix) + 10) / 2.54,
    height = (nrow(tf_cluster_matrix) + 10) / 2.54)

combined_heatmaps <- HeatmapList(heatmap_TFs + heatmap_ppa)
print(combined_heatmaps)
invisible(dev.off())

# Loading the gene reference data from genecode files
prom_th_up = if (nzchar(args[["prom_th_up"]])) {
  as.integer(args[["prom_th_up"]])
} else {
  message("Using default prom_th_up value: ", defaults$prom_th_up)
  defaults$prom_th_up
} 

prom_th_down = if (nzchar(args[["prom_th_down"]])) {
  as.integer(args[["prom_th_down"]])
} else {
  message("Using default prom_th_down value: ", defaults$prom_th_down)
  defaults$prom_th_down
}

gene_annot_dir = args[["gencode_dir"]]
gencode_transcripts  <- import(gene_annot_dir, format = "gff3", feature.type = "transcript")
gene_transcript_data <- gencode_transcripts[gencode_transcripts$transcript_type == "protein_coding"]
gene_id_lists <- mcols(gene_transcript_data)$gene_name 

pos_transcripts <- gene_transcript_data[strand(gene_transcript_data) == "+"]
neg_transcripts <- gene_transcript_data[strand(gene_transcript_data) == "-"]

prom_pos <- GRanges(
  seqnames = seqnames(pos_transcripts),
  ranges = IRanges(
    start = start(pos_transcripts) - prom_th_up,
    end   = start(pos_transcripts) + prom_th_down
  ),
  strand = strand(pos_transcripts),
  gene_name = mcols(pos_transcripts)$gene_name
)

prom_neg <- GRanges(
  seqnames = seqnames(neg_transcripts),
  ranges = IRanges(
    start = end(neg_transcripts) - prom_th_down,
    end   = end(neg_transcripts) + prom_th_up
  ),
  strand = strand(neg_transcripts),
  gene_name = mcols(neg_transcripts)$gene_name
)

gene_tss_grg <- c(prom_pos, prom_neg)

#Save this granges object as an RDS
saveRDS(gene_tss_grg, paste0(output_file_main,"gene_tss_granges.rds"))

#Directly overlapping genes' promoters with risk-mediating peaks
peak_ppa_frame_filt = peak_ppa_frame
peak_ppa_frame_filt$sumPPA = NULL
peak_ppa_frame_filt = peak_ppa_frame_filt %>% distinct()

rmp_granges = StringToGRanges(peak_ppa_frame_filt$region)

rmp_promoter_overlap = findOverlaps(rmp_granges, gene_tss_grg)

#Creating a dataframe to record these results
if (length(rmp_promoter_overlap) != 0) {
  direct_overlap_df = as.data.frame(matrix(0, nrow = length(rmp_promoter_overlap),
                                            ncol = 4))
  colnames(direct_overlap_df) = c("RMP","Promoter",
                                  "Gene","Cell_Type")
  
  direct_overlap_df$RMP = peak_ppa_frame_filt$region[queryHits(rmp_promoter_overlap)]
  direct_overlap_df$Promoter = GRangesToString(gene_tss_grg[subjectHits(rmp_promoter_overlap)])
  direct_overlap_df$Gene = gene_tss_grg$gene_name[subjectHits(rmp_promoter_overlap)]
  direct_overlap_df$Cell_Type = peak_ppa_frame_filt$cell[queryHits(rmp_promoter_overlap)]
  direct_overlap_df = direct_overlap_df %>%
    separate_rows(Cell_Type, sep = ',')
  
  #Output this as a spreadsheet
  direct_overlap_dir = paste0(output_file_main, "direct_rmp_gene_overlaps.txt")
  write.table(direct_overlap_df, file = direct_overlap_dir ,row.names = FALSE, quote = FALSE, sep = "\t")
} else {
  message('No genes found by direct overlap of RMPs with promoter peaks.')
}

#Getting the list of all Cicero files in the working directory and
#loading them and creating the list of cell types based on the file names
cicero_file_list = list.files(path = output_file_main, pattern = "filtered_coaccessible_sites4s.txt",
                              full.names = FALSE)
cell_type_list = sub(".filtered_coaccessible_sites4s.txt","",cicero_file_list)
cicero_file_list = paste0(output_file_main, cicero_file_list)

#Looping through the cell types, reading the Cicero file of that cell type
#and linking them to the Effect SNPs and gene promoters
coaccess_gene_cell_final = c()
for (i in c(1:length(cell_type_list))){
  #Getting the current cell type name and its Cicero file
  cell_temp = cell_type_list[i]
  cell_cicero_data = read.table(file = cicero_file_list[i],
                                sep = "\t", header = TRUE)
                                
  #Finding the peaks in the first column of the cell type data 
  #which overlaps the Effect SNPs and only keeping those peaks
  cell_peak1_grg = StringToGRanges(gsub("_","-",cell_cicero_data$Peak1))
  peak1_snp_overlap = findOverlaps(cell_peak1_grg, eff_snp_grg)
  cell_cicero_data = cell_cicero_data[unique(queryHits(peak1_snp_overlap)),]

  #Now mapping the peak2 of the cell type filtered Cicero data to the
  #promoter regions of genes
  cell_peak2_grg = StringToGRanges(gsub("_","-",cell_cicero_data$Peak2))
  peak2_promoter_overlap = findOverlaps(cell_peak2_grg, gene_tss_grg)
  
  #Creating a new dataframe to hold the peak1, peak2, coaccess,
  #promoter regions, and gene names for the list of peak2 ones
  #that overlap a promoter regions
  peak_gene_frame_temp = as.data.frame(matrix(nrow = length(peak2_promoter_overlap),
                                              ncol = 6))
  colnames(peak_gene_frame_temp) = c("Peak1","Peak2","coaccess","Promoter",
                                     "Gene","Cell_Type")
  
  if(nrow(peak_gene_frame_temp) == 0) { # check if no results for this cell type
    next
  }
  peak_gene_frame_temp$Peak1 = cell_cicero_data$Peak1[queryHits(peak2_promoter_overlap)]
  peak_gene_frame_temp$Peak2 = cell_cicero_data$Peak2[queryHits(peak2_promoter_overlap)]
  peak_gene_frame_temp$coaccess = cell_cicero_data$coaccess[queryHits(peak2_promoter_overlap)]
  peak_gene_frame_temp$Gene = gene_tss_grg$gene_name[subjectHits(peak2_promoter_overlap)]
  peak_gene_frame_temp$Promoter = GRangesToString(gene_tss_grg[subjectHits(peak2_promoter_overlap)])
  peak_gene_frame_temp$Cell_Type = cell_temp
  
  coaccess_gene_cell_final[[cell_temp]] = peak_gene_frame_temp
}

#The dataframe having the list of SNP-affected links matching 
#Promoter regions of genes for each cell type
#Checking to see if there were any results for cicero
if (is.null(coaccess_gene_cell_final)) {
  message('No genes found by overlap of effect-SNPs with putative enhancers using Cicero.')
} else {
  coaccess_gene_cell_final = do.call(rbind,coaccess_gene_cell_final)
  coaccess_gene_cell_final$Peak1 = gsub("_","-",coaccess_gene_cell_final$Peak1)
  coaccess_gene_cell_final$Peak2 = gsub("_","-",coaccess_gene_cell_final$Peak2)
  
  #Saving the final table
  cic_peak_interact_dir = paste0(output_file_main, "cic_peak_interact_gene.txt")
  write.table(coaccess_gene_cell_final, file = cic_peak_interact_dir, sep= "\t", row.names = FALSE, quote = FALSE)
}

#Filtering the dataframe to only include gene-cell data
gene_cell_final1 = coaccess_gene_cell_final %>%
  select(Cell_Type, Gene)

gene_cell_final2 = direct_overlap_df %>%
  select(Cell_Type, Gene)

gene_cell_final <- rbind(gene_cell_final1, gene_cell_final2)

#Creating a cell by gene matrix 
gene_names = unique(gene_cell_final$Gene)
cell_names = unique(gene_cell_final$Cell_Type)

gene_cell_matrix = matrix(0, nrow = length(cell_names), ncol = length(gene_names))
row.names(gene_cell_matrix) = cell_names
colnames(gene_cell_matrix) = gene_names

for (i in 1:nrow(gene_cell_final)) {
  row <- gene_cell_final$Cell_Type[i]
  col <- gene_cell_final$Gene[i]
  gene_cell_matrix[row, col] <- 1
}

#Saving the final table and matrix and a heatmap
cic_peak_interact_dir = paste0(output_file_main, "cic_peak_interact_gene.txt")
write.table(coaccess_gene_cell_final, file = cic_peak_interact_dir, sep="\t", row.names = FALSE, quote= FALSE)

cell_gene_out = paste0(output_file_main, "cell_gene_matrix.txt")
write.table(gene_cell_matrix, file = cell_gene_out, row.names = T, quote = F, sep = '\t')

if (length(rmp_promoter_overlap) != 0) {
  promoter_cell_gene = direct_overlap_df
  promoter_cell_gene[["Peak_gene"]] = paste0(promoter_cell_gene$RMP,"; ",
                                             promoter_cell_gene$Gene)
  
  prom_matrix = promoter_cell_gene[,c("Peak_gene","Cell_Type")]
  prom_matrix = unique(prom_matrix)
  prom_matrix <- prom_matrix %>%
    mutate(value = 1) %>%
    spread(key = Cell_Type, value = value, fill = 0)
  prom_matrix = tibble::column_to_rownames(prom_matrix, var = "Peak_gene")
  prom_matrix = as.matrix(prom_matrix)
  
  prom_matrix_file = paste0(output_file_main,"Gene_promoter_matrix.txt")
  write.table(prom_matrix, prom_matrix_file, sep = "\t", row.names = TRUE, quote = FALSE)
}

if (!is.null(coaccess_gene_cell_final)) {
  enh_cell_gene = coaccess_gene_cell_final
  enh_cell_gene[["Peak_gene"]] = paste0(enh_cell_gene$Peak1,"; ",
                                        enh_cell_gene$Gene)     
  
  enh_matrix = enh_cell_gene[,c("Peak_gene","Cell_Type")]
  enh_matrix = unique(enh_matrix)
  enh_matrix <- enh_matrix %>%
    mutate(value = 1) %>%
    spread(key = Cell_Type, value = value, fill = 0)
  enh_matrix = tibble::column_to_rownames(enh_matrix, var = "Peak_gene")
  enh_matrix = as.matrix(enh_matrix)
  
  enh_matrix_file = paste0(output_file_main,"Gene_enhancer_matrix.txt")
  write.table(enh_matrix, enh_matrix_file, sep = "\t", row.names = TRUE, quote = FALSE)
}
