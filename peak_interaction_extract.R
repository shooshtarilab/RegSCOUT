library(Repitools)
library(Seurat)
library(Signac)
library(monocle3)
library(cicero)
library(xlsx)
library(stringr)
library(circlize)
library(ComplexHeatmap)
library(R.utils)
library(ape)


#Setting command line arguments
args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Defining default parameter values
defaults <- list(
  coaccess_th = 0.05,
  cic_genomic_window = 2000000,
  peak_th = 0.1
)

#Getting the genome built of the scATAC-seq data
genome_built = args[["genome_built"]]

#Loading the reference genome object based on the genome built of the data 
if (genome_built == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
  genome = "BSgenome.Hsapiens.UCSC.hg19"  
  hg_19_ref = BSgenome.Hsapiens.UCSC.hg19@seqinfo
  hg_19_dat = as.data.frame(matrix(NA, nrow = length(hg_19_ref@seqnames), ncol = 2))
  hg_19_dat$V1 = hg_19_ref@seqnames
  hg_19_dat$V2 = hg_19_ref@seqlengths
}else if(genome_built == "hg38"){
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome = "BSgenome.Hsapiens.UCSC.hg38" 
  hg_38_ref = BSgenome.Hsapiens.UCSC.hg38@seqinfo
  hg_38_dat = as.data.frame(matrix(NA, nrow = length(hg_38_ref@seqnames), ncol = 2))
  hg_38_dat$V1 = hg_38_ref@seqnames
  hg_38_dat$V2 = hg_38_ref@seqlengths
}else if(genome_built == "mm9"){
  library(BSgenome.Mmusculus.UCSC.mm9)
  genome = "BSgenome.Mmusculus.UCSC.mm9" 
  mm_9_ref = BSgenome.Mmusculus.UCSC.mm9@seqinfo
  mm_9_dat = as.data.frame(matrix(NA, nrow = length(mm_9_ref@seqnames), ncol = 2))
  mm_9_dat$V1 = mm_9_ref@seqnames
  mm_9_dat$V2 = mm_9_ref@seqlengths
}else if(genome_built == "mm10"){
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome = "BSgenome.Mmusculus.UCSC.mm10"
  mm_10_ref = BSgenome.Mmusculus.UCSC.mm10@seqinfo
  mm_10_dat = as.data.frame(matrix(NA, nrow = length(mm_10_ref@seqnames), ncol = 2))
  mm_10_dat$V1 = mm_10_ref@seqnames
  mm_10_dat$V2 = mm_10_ref@seqlengths
}

#Loading the scATAC-seq object 
pbmc_dir = args[["seurat_obj"]]
pbmc_name = load(pbmc_dir) 
pbmc = get(pbmc_name)
rm(pbmc_name)
print("Seurat object loaded")

#Getting the output directory
output_file_main = args[["output_dir"]]

#Loading the cell information from the metadata of the scATAC-seq object
pbmc_meta = pbmc@meta.data
cell_info = rownames(pbmc_meta)
cellinfo = as.data.frame(cell_info)
names(cellinfo) = "cells"

#Loading the peak information from the peaks assay of the scATAC-seq object
pbmc_assays = pbmc@assays
peaks_assays = pbmc_assays$peaks
peaks_granges = peaks_assays@ranges
peaks_datatframe = annoGR2DF(peaks_granges)
peakinfo = subset(peaks_datatframe, select = -c(width))
names(peakinfo) = c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="-")
row.names(peakinfo) <- peakinfo$site_name

#Extracting the counts matrix and binarizing it
indata = peaks_assays@counts
indata@x[indata@x > 0] <- 1
row.names(indata) = row.names(peakinfo)
colnames(indata) = cellinfo$cells

#Preparing the counts matrix, peaks and cells dataframes for running the Cicero
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
row.names(fd) = row.names(indata)
row.names(pd) = colnames(indata)


#Running the Cicero in each cell type separately
print("Cell type specific loop started!")

coaccess_th = if (!is.null(args[["coaccess_th"]])) as.numeric(args[["coaccess_th"]]) else defaults$coaccess_th
cic_genomic_window = if (!is.null(args[["cic_genomic_window"]])) as.integer(args[["cic_genomic_window"]]) else defaults$cic_genomic_window
th = if (!is.null(args[["peak_th"]])) as.numeric(args[["peak_th"]]) else defaults$peak_th

for (cell in levels(pbmc)){
  #Extracting the peaks in the current cell type and filtering the matrix
  next_stat <- FALSE
  
  #Extracting the cells belonging to the current cell type from the counts matrix
  cell_idx = Idents(pbmc) == cell
  indata_cell = indata[,cell_idx]
  
  #Filtering the non relevant peaks in the data of the current cell type
  indata_temp = indata_cell[rowSums(indata_cell)>th*ncol(indata_cell),]
  
  
  print_message= paste0("number of rows=", nrow(indata_temp),
                        "number of columns=", ncol(indata_temp),
                        "maximum val = ", max(indata_temp))
  print(print_message)
  
  if(nrow(indata_temp) == 0){
    next
  }
  
  
  #Generating the cell dataset object for the cell type specific 
  #part of the counts matrix and running the Cicero
  
  tryCatch(
    expr = {
      input_cds <-  suppressWarnings(new_cell_data_set(indata_temp))
      
      set.seed(2017)
      input_cds <- detect_genes(input_cds)
      input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 
      input_cds = input_cds[,Matrix::colSums(exprs(input_cds)) != 0]
      input_cds <- estimate_size_factors(input_cds)
      input_cds <- preprocess_cds(input_cds, method = "PCA")
      input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP',
                                    preprocess_method = "PCA")
      
      umap_coords <- reducedDims(input_cds)$UMAP
      cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords,
                                    k = 30)
      if (genome_built == "hg19"){
        conns <- run_cicero(cicero_cds, genomic_coords = hg_19_dat, sample_num = 10000,
                            window = cic_genomic_window)
        
      }else if(genome_built == "hg38"){
        conns <- run_cicero(cicero_cds, genomic_coords = hg_38_dat, sample_num = 10000,
                            window = cic_genomic_window)
        
      }else if(genome_built == "mm9"){
        conns <- run_cicero(cicero_cds, genomic_coords = mm_9_dat, sample_num = 10000,
                            window = cic_genomic_window)
        
      }else if(genome_built == "mm10"){
        conns <- run_cicero(cicero_cds, genomic_coords = mm_10_dat, sample_num = 10000,
                            window = cic_genomic_window)
      }},
    error = function(e){
      print("Error happened!")
      next_stat <<- TRUE
    })
  print("Next stat value is:")
  print(next_stat)
  if(next_stat){
    next_stat <- FALSE
    print("Signal not found. Moving to the next cell")
    next
  }
  
  #Filtering the coaccessibility frame to remove rows with NA values
  #and the ones with coaccessibility lower than coaccess_th
  rem_index = !is.na(conns$coaccess)
  conns = conns[rem_index,]
  head(conns)

  conns$Peak1 = gsub("-","_",conns$Peak1)
  conns$Peak2 = gsub("-","_",conns$Peak2)
  conns_filt = conns[conns$coaccess>coaccess_th,]
  conns_filt = conns_filt[!is.na(conns_filt$Peak1),]
  
  if(nrow(conns_filt) == 0){
    next
  }
  
  #Saving the coaccessibitliy data for the current cell type
  conns_file = paste0(output_file_main,cell,".filtered_coaccessible_sites4s.txt")
  write.table(conns_filt, file = conns_file, sep = "\t", 
              col.names = TRUE, row.names = FALSE, quote = FALSE)
  print(cell)
}

print("Cicero is done and cell type specific coaccessibility tables are saved!")
