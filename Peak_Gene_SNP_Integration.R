library(GenomicRanges)
library(ggplot2)
library(stringr)
library(readxl)
library(dplyr)
library(tidyr)
library(xlsx)
library(ComplexHeatmap)
library(R.utils)
library(Seurat)
library(Signac)
library(ape)
library(tibble)
library(circlize)


#Setting command line arguments
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Getting the working directory
output_file_main = args[["output_dir"]]
#output_file_main = "~/rprojects_whole/Final_new_pip/results/MS/10x/"
print("Step passed")
#Getting the cell by peak table
cell_peak_file = paste0(output_file_main,"cell_peak.xlsx")
cell_peak = read_xlsx(path = cell_peak_file)

#Getting the file of effect SNPs and loading them
#eff_snp_file = args[["EffSNP_file"]]
#eff_snp_file = "~/rprojects_whole/ms_pbmc/gwas_data/bmi_chip/Ci_effect_SNPs.txt"
eff_snp_file = paste0(output_file_main,"Ci_effect_SNPs.txt")
eff_snp = read.table(file = eff_snp_file, header = TRUE)

#Finding the sum of PPA of peaks

peak_list = paste(cell_peak$chr,cell_peak$start,cell_peak$end,
                  sep = "-")
peak_granges = StringToGRanges(peak_list)
eff_snp_list = paste(eff_snp$CHR,eff_snp$Pos,(eff_snp$Pos+1),
                     sep = "-")
eff_snp_granges = StringToGRanges(eff_snp_list)

peak_eff_overlaps = findOverlaps(peak_granges,eff_snp_granges)

peak_ppa_frame = as.data.frame(matrix(0,nrow = length(peak_eff_overlaps),
                                      ncol = 3))
colnames(peak_ppa_frame) = c("region","PPA","cell")

peak_ppa_frame$region = peak_list[queryHits(peak_eff_overlaps)]
peak_ppa_frame$cell = cell_peak$`cell sub-types`[queryHits(peak_eff_overlaps)]
peak_ppa_frame$PPA = eff_snp$ppa[subjectHits(peak_eff_overlaps)]

peak_ppa_frame = peak_ppa_frame %>%
  group_by(region, cell) %>%
  mutate(sumPPA = sum(PPA)) %>%
  ungroup()

peak_ppa_frame = peak_ppa_frame %>% select(-PPA)
peak_ppa_frame = unique(peak_ppa_frame)

peak_ppa_frame_dir = paste0(output_file_main,"risk_regions_ppa.xlsx")
write.xlsx(as.data.frame(peak_ppa_frame), file = peak_ppa_frame_dir, col.names = TRUE,
           row.names = FALSE)

peak_cluster_matrix = peak_ppa_frame[,c("region","cell")]
peak_cluster_matrix <- peak_cluster_matrix %>%
  separate_rows(cell, sep = ",")
peak_cluster_matrix <- peak_cluster_matrix %>%
  mutate(value = 1) %>%
  spread(key = cell, value = value, fill = 0)
peak_cluster_matrix = tibble::column_to_rownames(peak_cluster_matrix, var = "region")
peak_cluster_matrix = as.matrix(peak_cluster_matrix)

peak_cluster_matrix_file = paste0(output_file_main,"peak_cluster_matrix.txt")
write.table(peak_cluster_matrix, peak_cluster_matrix_file, sep = "\t", row.names = TRUE)

f1 = colorRamp2(seq(0, 1, length = 2), c("#EEEEEE", "red"))


heatmap_peaks <- Heatmap(peak_cluster_matrix, name = "Percentage accessibility", col = f1, 
                         column_title = "Immune Cell Types", row_title = "Risk-mediating Peaks",
                         row_names_gp = grid::gpar(fontsize = 4),
                         column_names_gp = grid::gpar(fontsize = 15),
                         rect_gp = gpar(col= "#84878a"),
                         heatmap_legend_param = list(title = "Accessibility", at = c(0, 1), 
                                                     labels = c("0", "1"), 
                                                     color_bar = "vertical", 
                                                     legend_width = unit(16, "cm")))


# guideline bar based on sumPPA
risk_regions = peak_ppa_frame[,c("region", "sumPPA")]
rownames(risk_regions) = NULL
risk_regions = tibble::column_to_rownames(risk_regions, var = "region")
risk_regions = as.matrix(risk_regions)

# sorting the rownames of the risk_regions matrix to match the rownames of the peak_cluster_matrix
roworder <- order(match(rownames(risk_regions), rownames(peak_cluster_matrix)))
risk_regions <- risk_regions[roworder, , drop = FALSE]
#rownames(risk_regions) == rownames(peak_cluster_matrix)

f2 = colorRamp2(seq(0, max(risk_regions[,1], na.rm = TRUE), length = 2), c("#EEEEEE","blue"))

# create a separate heatmap for the ppa
heatmap_ppa <- Heatmap(risk_regions, name = "SumPPA", col = f2, 
                       heatmap_legend_param = list(title = "SumPPA", 
                                                   at = c(0, max(risk_regions, na.rm = TRUE)), 
                                                   labels = c("0", round(max(risk_regions, na.rm = TRUE))), 
                                                   color_bar = "vertical", 
                                                   legend_width = unit(16, "cm")),
                       row_names_gp = grid::gpar(fontsize = 4),
                       column_names_gp = grid::gpar(fontsize = 15),
                       # Add cell_fun argument to display values inside the heatmap cells
                       layer_fun = function(j, i, x, y, width, height, fill) {
                         grid::grid.text(sprintf("%.3f", risk_regions[i, j]), x, y, gp = grid::gpar(fontsize = 4))
                       })

# combine the two heatmaps and save it 
output_peak_file = paste0(output_file_main,"cell_peak.png")
png(output_peak_file, width = 2400, height = 8000, res = 300)
combined_heatmaps <- HeatmapList(heatmap_peaks + heatmap_ppa)
print(combined_heatmaps)
dev.off()

#Changing the cell by peak matrix so that each line has only one cell type
cell_peak = cell_peak %>% 
  separate_rows('cell sub-types', sep=",")


#Getting the peaks and Effect SNPs as GRanges and overlapping them to find
#SNP-affected peaks and mapping TFs to those peaks 
cell_peak_grg = StringToGRanges(paste(cell_peak$chr,cell_peak$start,
                                      cell_peak$end,sep = "-"))
eff_snp_grg = StringToGRanges(paste(eff_snp$CHR,(eff_snp$Pos-1),
                                    eff_snp$Pos,sep = "-"))


peak_snp_overlap = findOverlaps(cell_peak_grg, eff_snp_grg)


cell_peak_filt = as.data.frame(matrix(0,nrow = length(peak_snp_overlap),
                                      ncol=(ncol(cell_peak)+4)))

colnames(cell_peak_filt) = c(colnames(cell_peak),"SNP","TF","Locus","log_lik_ratio")

cell_peak_filt[,colnames(cell_peak)] = cell_peak[queryHits(peak_snp_overlap),]

cell_peak_filt$SNP = eff_snp$SNP[subjectHits(peak_snp_overlap)]
cell_peak_filt$TF = eff_snp$TF[subjectHits(peak_snp_overlap)]
cell_peak_filt$Locus = eff_snp$Locus[subjectHits(peak_snp_overlap)]
cell_peak_filt$log_lik_ratio = abs(eff_snp$log_like_ratio[subjectHits(peak_snp_overlap)])

cell_tf_snp = as.data.frame(matrix(0, nrow = nrow(cell_peak_filt),
                                   ncol = 4))
colnames(cell_tf_snp) = c("region","cell","TFSNP","log_lik_ratio")
cell_tf_snp$region = paste(cell_peak_filt$chr,
                           cell_peak_filt$start,
                           cell_peak_filt$end,
                           sep = "-")
cell_tf_snp$cell = cell_peak_filt$`cell sub-types`
cell_tf_snp$TFSNP = paste0(cell_peak_filt$TF,"-",cell_peak_filt$SNP)
cell_tf_snp$log_lik_ratio = cell_peak_filt$log_lik_ratio

peak_ratio_frame_dir = paste0(output_file_main,"risk_regions_ratio.xlsx")
write.xlsx(as.data.frame(cell_tf_snp), file = peak_ratio_frame_dir, col.names = TRUE,
           row.names = FALSE)

tf_cluster_matrix = cell_tf_snp[,c("TFSNP","cell")]
#tf_cluster_matrix <- tf_cluster_matrix %>%
#  separate_rows(cell, sep = ",")
tf_cluster_matrix <- tf_cluster_matrix %>%
  mutate(value = 1) %>%
  spread(key = cell, value = value, fill = 0)
tf_cluster_matrix = tibble::column_to_rownames(tf_cluster_matrix, var = "TFSNP")
tf_cluster_matrix = as.matrix(tf_cluster_matrix)

TF_cluster_matrix_file = paste0(output_file_main,"TF_cluster_matrix.txt")
write.table(tf_cluster_matrix, TF_cluster_matrix_file, sep = "\t", row.names = TRUE)


f1 = colorRamp2(seq(0, 1, length = 2), c("#EEEEEE", "red"))


heatmap_TFs <- Heatmap(tf_cluster_matrix, name = "Percentage accessibility", col = f1, 
                         column_title = "Immune Cell Types", row_title = "Risk-mediating TFs",
                         row_names_gp = grid::gpar(fontsize = 4),
                         column_names_gp = grid::gpar(fontsize = 15),
                         rect_gp = gpar(col= "#84878a"),
                         heatmap_legend_param = list(title = "Accessibility", at = c(0, 1), 
                                                     labels = c("0", "1"), 
                                                     color_bar = "vertical", 
                                                     legend_width = unit(16, "cm")))


# guideline bar based on sumPPA
risk_tfs = cell_tf_snp[,c("TFSNP", "log_lik_ratio")]
risk_tfs = unique(risk_tfs)
rownames(risk_tfs) = NULL
risk_tfs = tibble::column_to_rownames(risk_tfs, var = "TFSNP")
risk_tfs = as.matrix(risk_tfs)

# sorting the rownames of the risk_regions matrix to match the rownames of the peak_cluster_matrix
roworder <- order(match(rownames(risk_tfs), rownames(tf_cluster_matrix)))
risk_tfs <- risk_tfs[roworder, , drop = FALSE]
#rownames(risk_regions) == rownames(peak_cluster_matrix)

f2 = colorRamp2(seq(0, max(risk_tfs[,1], na.rm = TRUE), length = 2), c("#EEEEEE","blue"))

# create a separate heatmap for the ppa
heatmap_ppa <- Heatmap(risk_tfs, name = "Log_lik_ratio", col = f2, 
                       heatmap_legend_param = list(title = "Log_lik_ratio", 
                                                   at = c(0, max(risk_tfs, na.rm = TRUE)), 
                                                   labels = c("0", round(max(risk_tfs, na.rm = TRUE))), 
                                                   color_bar = "vertical", 
                                                   legend_width = unit(16, "cm")),
                       row_names_gp = grid::gpar(fontsize = 4),
                       column_names_gp = grid::gpar(fontsize = 15),
                       # Add cell_fun argument to display values inside the heatmap cells
                       layer_fun = function(j, i, x, y, width, height, fill) {
                         grid::grid.text(sprintf("%.3f", risk_tfs[i, j]), x, y, gp = grid::gpar(fontsize = 4))
                       })

# combine the two heatmaps and save it 
output_tf_file = paste0(output_file_main,"cell_tf.png")
png(output_tf_file, width = 2400, height = 8000, res = 300)
combined_heatmaps <- HeatmapList(heatmap_TFs + heatmap_ppa)
print(combined_heatmaps)
dev.off()


#In the final dataframe, the affected peaks, SNPs affecting them,
#the loci of SNPs, and TFs are provided

# cell_peak_filt[["locs-reg"]] = paste0("Loc ",cell_peak_filt$Locus,":",
#                                       paste(cell_peak_filt$chr,
#                                             cell_peak_filt$start,
#                                             cell_peak_filt$end,
#                                             sep = "-"))
# 
# 
# 
# #Creating cell by TF and cell by peak-loci matrices
# cell_tf_mat = as.matrix(table(cell_peak_filt$`cell sub-types`, cell_peak_filt$TF))
# cell_tf_mat[which(cell_tf_mat > 0)] = 1
# 
# cell_peak_mat = as.matrix(table(cell_peak_filt$`cell sub-types`, cell_peak_filt$`locs-reg`))
# cell_peak_mat[which(cell_peak_mat > 0)] = 1
# 
# #Saving the matrices and generating heatmaps of results
# cell_peak_out = paste0(output_file_main, "filt_cell_peak.csv")
# write.csv(cell_peak_filt, file = cell_peak_out, row.names = FALSE,
#           col.names = TRUE, quote = FALSE)
# 
# cell_tf_out = paste0(output_file_main, "cell_tf_mat.csv")
# write.csv(cell_tf_mat, file = cell_tf_out, row.names = TRUE,
#           col.names = TRUE, quote = FALSE)
# 
# cell_loc_out = paste0(output_file_main, "cell_reg.csv")
# write.csv(cell_peak_mat, file = cell_loc_out, row.names = TRUE,
#           col.names = TRUE, quote = FALSE)
# 
# 
# f1 = c("0" = "white", "1" = "blue")
# 
# output_file = paste0(output_file_main, "cell_tf.png")
# file.remove(output_file)
# png(output_file,width = 6400, height = 2800, res = 300)
# Heatmap(cell_tf_mat, name = "TF presence", col = f1, 
#         column_title = "TF-cell type plot",
#         row_names_gp = grid::gpar(fontsize = 10),
#         column_names_gp = grid::gpar(fontsize = 10),
#         rect_gp = gpar(col= "#84878a"))
# 
# dev.off()
# 
# output_file = paste0(output_file_main, "cell_peak_locus.png")
# file.remove(output_file)
# png(output_file,width = 8400, height = 3800, res = 300)
# Heatmap(cell_peak_mat, name = "Peak presence", col = f1, 
#         column_title = "Peak-locus-cell type plot",
#         row_names_gp = grid::gpar(fontsize = 10),
#         column_names_gp = grid::gpar(fontsize = 10),
#         rect_gp = gpar(col= "#84878a"))
# 
# dev.off()

print("Affected peaks and TFs extraction finished!")

#Loading the gene reference data from genecode files
prom_th_up = as.numeric(args[["prom_th_up"]])
prom_th_down = as.numeric(args[["prom_th_down"]])
gene_annot_dir = args[["genecode_dir"]]
#gene_annot_dir = "~/rprojects_whole/ATAC_seq_RA/gwas_data/naderH_new/nadeR/files/gencode.v43lift37.annotation.gff3"
gene_annot = read.gff(gene_annot_dir, na.strings = c(".", "?"), GFF3 = TRUE)
gene_type_list = str_split(gene_annot$attributes, "gene_type=", simplify = TRUE)
gene_type_list = str_split(gene_type_list[,2], ";", simplify = TRUE)[,1]

gene_coding_index = gene_type_list == "protein_coding"
gene_annot = gene_annot[gene_coding_index,]
gene_transcript_data = gene_annot[gene_annot$type == "gene",]
gene_id_list = gene_transcript_data$attributes
gene_id_list = str_split(gene_id_list, "gene_name=", simplify = TRUE)
gene_id_list = str_split(gene_id_list[,2], ";", simplify = TRUE)[,1]

gene_transcript_data[["gene_name"]] = gene_id_list
pos_strand_index = gene_transcript_data$strand == "+"
neg_strand_index = gene_transcript_data$strand == "-"

gene_transcript_data[["TSS"]] = rep(0, nrow(gene_transcript_data))
gene_transcript_data$TSS[pos_strand_index] = gene_transcript_data$start[pos_strand_index]
gene_transcript_data$TSS[neg_strand_index] = gene_transcript_data$end[neg_strand_index]
gene_transcript_data[["length"]] = gene_transcript_data$end - gene_transcript_data$start

gene_id_list = str_split(gene_transcript_data$attributes, "gene_id=", simplify = TRUE)
gene_id_list = str_split(gene_id_list[,2], ";", simplify = TRUE)[,1]
gene_id_list = str_split(gene_id_list, "[.]", simplify = TRUE)[,1]
print(head(gene_id_list,10))

gene_transcript_data[["ens_id"]] = gene_id_list

gene_data_temp_pos = gene_transcript_data[gene_transcript_data$strand == "+",]

gene_tss_grg_pos = GRanges(ranges= IRanges(start = gene_data_temp_pos$TSS - prom_th_down,
                                       end = gene_data_temp_pos$TSS + prom_th_up),
                       seqnames = gene_data_temp_pos$seqid)

gene_tss_grg_pos$gene_name = gene_data_temp_pos$gene_name

gene_data_temp_neg = gene_transcript_data[gene_transcript_data$strand == "-",]

gene_tss_grg_neg = GRanges(ranges= IRanges(start = gene_data_temp_neg$TSS - prom_th_up,
                                           end = gene_data_temp_neg$TSS + prom_th_down),
                           seqnames = gene_data_temp_neg$seqid)

gene_tss_grg_neg$gene_name = gene_data_temp_neg$gene_name

gene_tss_grg = c(gene_tss_grg_pos, gene_tss_grg_neg)


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
  
  #Adding some lines having affected peaks repeated in both columns
  #with coaccessibility 1 to also capture genes with affected promotors
  cell_cicero_data_rep = as.data.frame(matrix(0, nrow = length(unique(cell_cicero_data$Peak1)),
                                              ncol = ncol(cell_cicero_data)))
  colnames(cell_cicero_data_rep) = colnames(cell_cicero_data)
  
  cell_cicero_data_rep$Peak1 = unique(cell_cicero_data$Peak1)
  cell_cicero_data_rep$Peak2 = unique(cell_cicero_data$Peak1)
  cell_cicero_data_rep$coaccess = 1
  
  cell_cicero_data = rbind(cell_cicero_data, cell_cicero_data_rep)
  #Now mapping the peak2 of the cell type filtered Cicero data to the
  #promotor regions of genes
  cell_peak2_grg = StringToGRanges(gsub("_","-",cell_cicero_data$Peak2))
  peak2_promotor_overlap = findOverlaps(cell_peak2_grg, gene_tss_grg)
  
  #Creating a new dataframe to hold the peak1, peak2, coaccess,
  #promotor regions, and gene names for the list of peak2 ones
  #that overlap a promotor regions
  peak_gene_frame_temp = as.data.frame(matrix(nrow = length(peak2_promotor_overlap),
                                              ncol = 6))
  colnames(peak_gene_frame_temp) = c("Peak1","Peak2","coaccess","Promotor",
                                     "Gene","Cell_Type")
  
  peak_gene_frame_temp$Peak1 = cell_cicero_data$Peak1[queryHits(peak2_promotor_overlap)]
  peak_gene_frame_temp$Peak2 = cell_cicero_data$Peak2[queryHits(peak2_promotor_overlap)]
  peak_gene_frame_temp$coaccess = cell_cicero_data$coaccess[queryHits(peak2_promotor_overlap)]
  peak_gene_frame_temp$Gene = gene_tss_grg$gene_name[subjectHits(peak2_promotor_overlap)]
  peak_gene_frame_temp$Promotor = GRangesToString(gene_tss_grg[subjectHits(peak2_promotor_overlap)])
  peak_gene_frame_temp$Cell_Type = cell_temp
  
  coaccess_gene_cell_final[[cell_temp]] = peak_gene_frame_temp
  
}

#The datafrae having the list of SNP-affected links matching 
#Promotor regions of genes for each cell type
coaccess_gene_cell_final = do.call(rbind,coaccess_gene_cell_final)

#Filtering the dataframe to only include gene-cell data
gene_cell_final = coaccess_gene_cell_final %>%
  select(Cell_Type, Gene)


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
peak_interact_gene_dir = paste0(output_file_main, "Aff_peak_interct_gene.csv")
write.csv(coaccess_gene_cell_final, file = peak_interact_gene_dir, row.names = FALSE,
          col.names = TRUE, quote = FALSE)

cell_gene_out = paste0(output_file_main, "cell_gene.csv")
write.csv(gene_cell_matrix, file = cell_gene_out, row.names = TRUE,
          col.names = TRUE, quote = FALSE)


f1 = c("0" = "white", "1" = "red")

output_file = paste0(output_file_main, "cell_gene.png")
file.remove(output_file)
png(output_file,width = 8400, height = 2800, res = 300)
Heatmap(gene_cell_matrix, name = "Gene presence", col = f1, 
        column_title = "Gene-cell type plot",
        row_names_gp = grid::gpar(fontsize = 16),
        column_names_gp = grid::gpar(fontsize = 6),
        rect_gp = gpar(col= "#84878a"))

dev.off()


promotor_cell_gene = coaccess_gene_cell_final[coaccess_gene_cell_final$Peak1 == coaccess_gene_cell_final$Peak2,]
promotor_cell_gene[["Peak_gene"]] = paste0(promotor_cell_gene$Peak1,"; ",
                                           promotor_cell_gene$Gene)

prom_matrix = promotor_cell_gene[,c("Peak_gene","Cell_Type")]
#tf_cluster_matrix <- tf_cluster_matrix %>%
#  separate_rows(cell, sep = ",")
prom_matrix <- prom_matrix %>%
  mutate(value = 1) %>%
  spread(key = Cell_Type, value = value, fill = 0)
prom_matrix = tibble::column_to_rownames(prom_matrix, var = "Peak_gene")
prom_matrix = as.matrix(prom_matrix)

prom_matrix_file = paste0(output_file_main,"Gene_promotor_matrix.txt")
write.table(prom_matrix, prom_matrix_file, sep = "\t", row.names = TRUE)

f1 = c("0" = "white", "1" = "red")

output_file = paste0(output_file_main, "cell_gene_promotor.png")
file.remove(output_file)
png(output_file,width = 8400, height = 2800, res = 300)
Heatmap(prom_matrix, name = "Gene presence", col = f1, 
        column_title = "Gene-cell type plot",
        row_names_gp = grid::gpar(fontsize = 6),
        column_names_gp = grid::gpar(fontsize = 16),
        rect_gp = gpar(col= "#84878a"))

dev.off()


enh_cell_gene = coaccess_gene_cell_final[coaccess_gene_cell_final$Peak1 != coaccess_gene_cell_final$Peak2,]
enh_cell_gene[["Peak_gene"]] = paste0(enh_cell_gene$Peak1,"; ",
                                           enh_cell_gene$Gene)

enh_matrix = enh_cell_gene[,c("Peak_gene","Cell_Type")]
enh_matrix = unique(enh_matrix)
#tf_cluster_matrix <- tf_cluster_matrix %>%
#  separate_rows(cell, sep = ",")
enh_matrix <- enh_matrix %>%
  mutate(value = 1) %>%
  spread(key = Cell_Type, value = value, fill = 0)
enh_matrix = tibble::column_to_rownames(enh_matrix, var = "Peak_gene")
enh_matrix = as.matrix(enh_matrix)

enh_matrix_file = paste0(output_file_main,"Gene_enhancer_matrix.txt")
write.table(enh_matrix, enh_matrix_file, sep = "\t", row.names = TRUE)

f1 = c("0" = "white", "1" = "red")

output_file = paste0(output_file_main, "cell_gene_enhancer.png")
file.remove(output_file)
png(output_file,width = 1400, height = 6800, res = 300)
Heatmap(enh_matrix, name = "Gene presence", col = f1, 
        column_title = "Gene-cell type plot",
        row_names_gp = grid::gpar(fontsize = 2),
        column_names_gp = grid::gpar(fontsize = 16),
        rect_gp = gpar(col= "#84878a"))

dev.off()







print("Pipeline finished!")
