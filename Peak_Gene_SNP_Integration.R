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
eff_snp_file = args[["EffSNP_file"]]
#eff_snp_file = "~/rprojects_whole/ms_pbmc/gwas_data/bmi_chip/Ci_effect_SNPs.txt"
eff_snp = read.table(file = eff_snp_file, header = TRUE)


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
                                      ncol=(ncol(cell_peak)+3)))

colnames(cell_peak_filt) = c(colnames(cell_peak),"SNP","TF","Locus")

cell_peak_filt[,colnames(cell_peak)] = cell_peak[queryHits(peak_snp_overlap),]

cell_peak_filt$SNP = eff_snp$SNP[subjectHits(peak_snp_overlap)]
cell_peak_filt$TF = eff_snp$TF[subjectHits(peak_snp_overlap)]
cell_peak_filt$Locus = eff_snp$Locus[subjectHits(peak_snp_overlap)]

#In the final dataframe, the affected peaks, SNPs affecting them,
#the loci of SNPs, and TFs are provided

cell_peak_filt[["locs-reg"]] = paste0("Loc ",cell_peak_filt$Locus,":",
                                      paste(cell_peak_filt$chr,
                                            cell_peak_filt$start,
                                            cell_peak_filt$end,
                                            sep = "-"))

#Creating cell by TF and cell by peak-loci matrices
cell_tf_mat = as.matrix(table(cell_peak_filt$`cell sub-types`, cell_peak_filt$TF))
cell_tf_mat[which(cell_tf_mat > 0)] = 1

cell_peak_mat = as.matrix(table(cell_peak_filt$`cell sub-types`, cell_peak_filt$`locs-reg`))
cell_peak_mat[which(cell_peak_mat > 0)] = 1

#Saving the matrices and generating heatmaps of results
cell_peak_out = paste0(output_file_main, "filt_cell_peak.csv")
write.csv(cell_peak_filt, file = cell_peak_out, row.names = FALSE,
          col.names = TRUE, quote = FALSE)

cell_tf_out = paste0(output_file_main, "cell_tf_mat.csv")
write.csv(cell_tf_mat, file = cell_tf_out, row.names = TRUE,
          col.names = TRUE, quote = FALSE)

cell_loc_out = paste0(output_file_main, "cell_reg.csv")
write.csv(cell_peak_mat, file = cell_loc_out, row.names = TRUE,
          col.names = TRUE, quote = FALSE)


f1 = c("0" = "white", "1" = "blue")

output_file = paste0(output_file_main, "cell_tf.png")
file.remove(output_file)
png(output_file,width = 6400, height = 2800, res = 300)
Heatmap(cell_tf_mat, name = "TF presence", col = f1, 
        column_title = "TF-cell type plot",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        rect_gp = gpar(col= "#84878a"))

dev.off()

output_file = paste0(output_file_main, "cell_peak_locus.png")
file.remove(output_file)
png(output_file,width = 8400, height = 3800, res = 300)
Heatmap(cell_peak_mat, name = "Peak presence", col = f1, 
        column_title = "Peak-locus-cell type plot",
        row_names_gp = grid::gpar(fontsize = 10),
        column_names_gp = grid::gpar(fontsize = 10),
        rect_gp = gpar(col= "#84878a"))

dev.off()

print("Affected peaks and TFs extraction finished!")

#Loading the gene reference data from genecode files
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

gene_tss_grg_pos = GRanges(ranges= IRanges(start = gene_data_temp_pos$TSS - 250,
                                       end = gene_data_temp_pos$TSS + 25),
                       seqnames = gene_data_temp_pos$seqid)

gene_tss_grg_pos$gene_name = gene_data_temp_pos$gene_name

gene_data_temp_neg = gene_transcript_data[gene_transcript_data$strand == "-",]

gene_tss_grg_neg = GRanges(ranges= IRanges(start = gene_data_temp_neg$TSS - 25,
                                           end = gene_data_temp_neg$TSS + 250),
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

print("Pipeline finished!")