suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(stringr))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(Signac))
suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(tools))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Function to read file based on extension
read_file <- function(file_path) {
  ext <- file_ext(file_path)
  
  # Decide based on extension
  if (ext == "csv") {
    df = read.csv(file_path, header = TRUE)
  } else if (ext == "tsv") {
    df = read.delim(file_path, header = TRUE)
  } else {
    df = read.table(file_path, header = TRUE)
  }
  return(df)
}

# read output directory
output_dir = args[["output_dir"]]

# defining default parameter values
defaults <- list(
  prom_th_up = 2000,
  prom_th_down = 2000
)

# read hic instructions directory
hic_instruct_dir <- args[["hic_instruct_dir"]]

# check if instructions spreadsheet exists
if (!file.exists(hic_instruct_dir)) {
  stop("Hi-C analysis requested but Hi-C instructions spreadsheet not found, please ensure path is correct/provided. Or if Hi-C analysis is not desired please do not use --hic_eqtl_analysis parameter.")
} 

# loading in user instructions
user_instruct <- read_file(hic_instruct_dir)
colnames(user_instruct) = tolower(colnames(user_instruct))

# load in gencode gene annotation and processing it
prom_th_up = if (nzchar(args[["prom_th_up"]])) {
  as.integer(args[["prom_th_up"]])
} else {
  defaults$prom_th_up
} 

prom_th_down = if (nzchar(args[["prom_th_down"]])) {
  as.integer(args[["prom_th_down"]])
} else {
  defaults$prom_th_down
}

gene_tss_grg = readRDS(paste0(output_dir, "gene_tss_granges.rds"))

# load in RMP information
rmp_df <- read.table(paste0(output_dir, 'risk_regions_ppa.txt'), header = TRUE)

# getting list of atac cell types requested and rmp cell types, seeing if rmps were not found in some atac cell types, removing them
rmp_cell_types <- unique(unlist(str_split(rmp_df$cell, ',')))
peak_cell_types <- unique(unlist(str_split(user_instruct$atac_cell_types, ',')))
missing_cell_types <- setdiff(peak_cell_types, rmp_cell_types)

if (length(missing_cell_types) > 0) { # removing missing cell types from user instructions
  message("Note: Hi-C analysis will not be conducted on cell types in which no SNPs co-localizing with risk-mediating peaks were identified: ", paste(missing_cell_types, collapse = ", "), ". Analysis not being conducted may also occur if a cell type name in the atac_cell_types column of the instructions spreadsheet does not match with the corresponding cell type name provided in the scATAC-seq object (i.e., B cell vs b_cell).")
}

# create separate chr, start, end columns for RMPs
rmp_df <- rmp_df %>%
  separate(region, into = c("chr", "start", "end"), sep = "-", convert = T)

# expand such that there is only one cell type per row for ease of filtering 
rmp_df <- rmp_df %>%
  separate_rows(cell, sep = ',') %>%
  distinct()

## defining functions for each type of Hi-C analysis (with/without gene information and single cell/bulk)
bulk_nogene_analysis <- function(hic_data, rmp_data, atac_ct, gencode_granges) {
  # filter hic dataset for certain columns, we assume all interactions in this Hi-C dataset are significant
  columns_to_keep <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2')
  hic_data <- hic_data[,columns_to_keep]
  
  # create two sets of granges, one for each set of genomic regions
  hic_granges1 <- GRanges(
    seqnames = hic_data$chr1,
    ranges = IRanges(start = hic_data$start1, end = hic_data$end1),
  )
  
  hic_granges2 <- GRanges(
    seqnames = hic_data$chr2,
    ranges = IRanges(start = hic_data$start2, end = hic_data$end2),
  )
  
  # overlap these granges with promoter granges as defined using gencode earlier
  hic_gene_overlap1 <- findOverlaps(hic_granges1, gencode_granges)
  hic_gene_overlap2 <- findOverlaps(hic_granges2, gencode_granges)
  
  # accounting for cases where there may be no overlaps b/w Hi-C data and gencode genes
  if (length(hic_gene_overlap1) == 0 & length(hic_gene_overlap2) == 0) {
    return('none')
  } else if (length(hic_gene_overlap1) == 0) {
    # working with just the hic_gene_overlap2 results
    gene_overlap_df2 <- data.frame(
      Chr1 = hic_data$chr1[queryHits(hic_gene_overlap2)],
      Start1 = hic_data$start1[queryHits(hic_gene_overlap2)],
      End1 = hic_data$end1[queryHits(hic_gene_overlap2)],
      promChr = hic_data$chr2[queryHits(hic_gene_overlap2)],
      promStart = hic_data$start2[queryHits(hic_gene_overlap2)],
      promEnd = hic_data$end2[queryHits(hic_gene_overlap2)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap2)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap2)])
    )
    
    # now create granges for other region that interacted with promoter region
    new_hic_granges2 <- GRanges(
      seqnames = gene_overlap_df2$Chr1,
      ranges = IRanges(start = gene_overlap_df2$Start1, end = gene_overlap_df2$End1),
    )
    
    # filter rmp_df based on atac_cell_types
    rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
    
    # create rmp granges
    rmp_granges <- GRanges(
      seqnames = rmp_df_filt$chr,
      ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
    )
    
    # find overlaps between these granges and RMPs
    hic_rmp_overlap2 <- findOverlaps(new_hic_granges2, rmp_granges)
    
    # accounting for case where no overlap b/w RMPs and the filtered Hi-C data
    if (length(hic_rmp_overlap2) == 0) {
      return('none')
    } else {
      # creating dataframes of results of this overlap
      rmp_gene_df2 <- data.frame(
        promChr = gene_overlap_df2$promChr[queryHits(hic_rmp_overlap2)],
        promStart = gene_overlap_df2$promStart[queryHits(hic_rmp_overlap2)],
        promEnd = gene_overlap_df2$promEnd[queryHits(hic_rmp_overlap2)],
        rmpChr = gene_overlap_df2$Chr1[queryHits(hic_rmp_overlap2)],
        rmpStart = gene_overlap_df2$Start1[queryHits(hic_rmp_overlap2)],
        rmpEnd = gene_overlap_df2$End1[queryHits(hic_rmp_overlap2)],
        Gene = gene_overlap_df2$Gene[queryHits(hic_rmp_overlap2)],
        Promoter = gene_overlap_df2$Promoter[queryHits(hic_rmp_overlap2)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap2)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap2)]
      )
      
      all_rmp_gene_df <- rmp_gene_df2
    }
  } else if (length(hic_gene_overlap2) == 0) {
    # create dataframes of results of this overlap
    gene_overlap_df1 <- data.frame(
      promChr = hic_data$chr1[queryHits(hic_gene_overlap1)],
      promStart = hic_data$start1[queryHits(hic_gene_overlap1)],
      promEnd = hic_data$end1[queryHits(hic_gene_overlap1)],
      Chr2 = hic_data$chr2[queryHits(hic_gene_overlap1)],
      Start2 = hic_data$start2[queryHits(hic_gene_overlap1)],
      End2 = hic_data$end2[queryHits(hic_gene_overlap1)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap1)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap1)])
    )
    
    # now create granges for other region that interacted with promoter region
    new_hic_granges1 <- GRanges(
      seqnames = gene_overlap_df1$Chr2,
      ranges = IRanges(start = gene_overlap_df1$Start2, end = gene_overlap_df1$End2),
    )
    
    # filter rmp_df based on atac_cell_types
    rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
    
    # create rmp granges
    rmp_granges <- GRanges(
      seqnames = rmp_df_filt$chr,
      ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
    )
    
    # find overlaps between these granges and RMPs
    hic_rmp_overlap1 <- findOverlaps(new_hic_granges1, rmp_granges)
    
    # accounting for case where no overlap b/w RMPs and the filtered Hi-C data
    if (length(hic_rmp_overlap1) == 0) {
      return('none')
    } else {
      # creating dataframes of results of this overlap
      rmp_gene_df1 <- data.frame(
        promChr = gene_overlap_df1$promChr[queryHits(hic_rmp_overlap1)],
        promStart = gene_overlap_df1$promStart[queryHits(hic_rmp_overlap1)],
        promEnd = gene_overlap_df1$promEnd[queryHits(hic_rmp_overlap1)],
        rmpChr = gene_overlap_df1$Chr2[queryHits(hic_rmp_overlap1)],
        rmpStart = gene_overlap_df1$Start2[queryHits(hic_rmp_overlap1)],
        rmpEnd = gene_overlap_df1$End2[queryHits(hic_rmp_overlap1)],
        Gene = gene_overlap_df1$Gene[queryHits(hic_rmp_overlap1)],
        Promoter = gene_overlap_df1$Promoter[queryHits(hic_rmp_overlap1)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap1)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap1)]
      )
      
      all_rmp_gene_df <- rmp_gene_df1 
    }
  } else {
    # create dataframes of results of this overlap
    gene_overlap_df1 <- data.frame(
      promChr = hic_data$chr1[queryHits(hic_gene_overlap1)],
      promStart = hic_data$start1[queryHits(hic_gene_overlap1)],
      promEnd = hic_data$end1[queryHits(hic_gene_overlap1)],
      Chr2 = hic_data$chr2[queryHits(hic_gene_overlap1)],
      Start2 = hic_data$start2[queryHits(hic_gene_overlap1)],
      End2 = hic_data$end2[queryHits(hic_gene_overlap1)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap1)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap1)])
    )
    
    gene_overlap_df2 <- data.frame(
      Chr1 = hic_data$chr1[queryHits(hic_gene_overlap2)],
      Start1 = hic_data$start1[queryHits(hic_gene_overlap2)],
      End1 = hic_data$end1[queryHits(hic_gene_overlap2)],
      promChr = hic_data$chr2[queryHits(hic_gene_overlap2)],
      promStart = hic_data$start2[queryHits(hic_gene_overlap2)],
      promEnd = hic_data$end2[queryHits(hic_gene_overlap2)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap2)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap2)])
    )
    
    # now create granges for other region that interacted with promoter region, keeping with the numbering convention of dfs
    new_hic_granges1 <- GRanges(
      seqnames = gene_overlap_df1$Chr2,
      ranges = IRanges(start = gene_overlap_df1$Start2, end = gene_overlap_df1$End2),
    )
    
    new_hic_granges2 <- GRanges(
      seqnames = gene_overlap_df2$Chr1,
      ranges = IRanges(start = gene_overlap_df2$Start1, end = gene_overlap_df2$End1),
    )
    
    # filter rmp_df based on atac_cell_types
    rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
    
    # create rmp granges
    rmp_granges <- GRanges(
      seqnames = rmp_df_filt$chr,
      ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
    )
    
    # find overlaps between these granges and RMPs
    hic_rmp_overlap1 <- findOverlaps(new_hic_granges1, rmp_granges)
    hic_rmp_overlap2 <- findOverlaps(new_hic_granges2, rmp_granges)
    
    # accounting for cases where there is no overlap between RMPs and the filtered Hi-C dataset
    if (length(hic_rmp_overlap1) == 0 & length(hic_rmp_overlap2) == 0) {
      return('none')
    } else if (length(hic_rmp_overlap1) == 0) {
      # creating dataframes of results of this overlap
      rmp_gene_df2 <- data.frame(
        promChr = gene_overlap_df2$promChr[queryHits(hic_rmp_overlap2)],
        promStart = gene_overlap_df2$promStart[queryHits(hic_rmp_overlap2)],
        promEnd = gene_overlap_df2$promEnd[queryHits(hic_rmp_overlap2)],
        rmpChr = gene_overlap_df2$Chr1[queryHits(hic_rmp_overlap2)],
        rmpStart = gene_overlap_df2$Start1[queryHits(hic_rmp_overlap2)],
        rmpEnd = gene_overlap_df2$End1[queryHits(hic_rmp_overlap2)],
        Gene = gene_overlap_df2$Gene[queryHits(hic_rmp_overlap2)],
        Promoter = gene_overlap_df2$Promoter[queryHits(hic_rmp_overlap2)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap2)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap2)]
      )
      
      all_rmp_gene_df <- rmp_gene_df2
    } else if (length(hic_rmp_overlap2) == 0) {
      # creating dataframes of results of this overlap
      rmp_gene_df1 <- data.frame(
        promChr = gene_overlap_df1$promChr[queryHits(hic_rmp_overlap1)],
        promStart = gene_overlap_df1$promStart[queryHits(hic_rmp_overlap1)],
        promEnd = gene_overlap_df1$promEnd[queryHits(hic_rmp_overlap1)],
        rmpChr = gene_overlap_df1$Chr2[queryHits(hic_rmp_overlap1)],
        rmpStart = gene_overlap_df1$Start2[queryHits(hic_rmp_overlap1)],
        rmpEnd = gene_overlap_df1$End2[queryHits(hic_rmp_overlap1)],
        Gene = gene_overlap_df1$Gene[queryHits(hic_rmp_overlap1)],
        Promoter = gene_overlap_df1$Promoter[queryHits(hic_rmp_overlap1)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap1)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap1)]
      )
      
      all_rmp_gene_df <- rmp_gene_df1
    } else {
      # creating dataframes of results of this overlap
      rmp_gene_df1 <- data.frame(
        promChr = gene_overlap_df1$promChr[queryHits(hic_rmp_overlap1)],
        promStart = gene_overlap_df1$promStart[queryHits(hic_rmp_overlap1)],
        promEnd = gene_overlap_df1$promEnd[queryHits(hic_rmp_overlap1)],
        rmpChr = gene_overlap_df1$Chr2[queryHits(hic_rmp_overlap1)],
        rmpStart = gene_overlap_df1$Start2[queryHits(hic_rmp_overlap1)],
        rmpEnd = gene_overlap_df1$End2[queryHits(hic_rmp_overlap1)],
        Gene = gene_overlap_df1$Gene[queryHits(hic_rmp_overlap1)],
        Promoter = gene_overlap_df1$Promoter[queryHits(hic_rmp_overlap1)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap1)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap1)]
      )
      
      rmp_gene_df2 <- data.frame(
        promChr = gene_overlap_df2$promChr[queryHits(hic_rmp_overlap2)],
        promStart = gene_overlap_df2$promStart[queryHits(hic_rmp_overlap2)],
        promEnd = gene_overlap_df2$promEnd[queryHits(hic_rmp_overlap2)],
        rmpChr = gene_overlap_df2$Chr1[queryHits(hic_rmp_overlap2)],
        rmpStart = gene_overlap_df2$Start1[queryHits(hic_rmp_overlap2)],
        rmpEnd = gene_overlap_df2$End1[queryHits(hic_rmp_overlap2)],
        Gene = gene_overlap_df2$Gene[queryHits(hic_rmp_overlap2)],
        Promoter = gene_overlap_df2$Promoter[queryHits(hic_rmp_overlap2)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap2)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap2)]
      )
      
      # bind these two dataframes tgt
      all_rmp_gene_df <- rbind(rmp_gene_df1, rmp_gene_df2)
    }
  }
  
  # remove all duplicate entries if present
  all_rmp_gene_df <- all_rmp_gene_df %>% distinct()
  
  # finalize this dataframe
  all_rmp_gene_df$promRegion <- paste(all_rmp_gene_df$promChr,
                                      all_rmp_gene_df$promStart, 
                                      all_rmp_gene_df$promEnd, sep = "-")
  
  all_rmp_gene_df$oeRegion <- paste(all_rmp_gene_df$rmpChr, 
                                    all_rmp_gene_df$rmpStart, 
                                    all_rmp_gene_df$rmpEnd, sep = "-")
  
  all_rmp_gene_df <- all_rmp_gene_df %>%
    select(-rmpChr, -rmpStart, -rmpEnd, -promChr, -promStart, -promEnd, -Promoter) %>%
    arrange(cell)
  
  colnames(all_rmp_gene_df)[colnames(all_rmp_gene_df) == 'Gene'] <- "gene"
  
  all_rmp_gene_df <- all_rmp_gene_df[,c('cell', 'gene', 'rmpRegion', 'promRegion', 'oeRegion')] %>% distinct()
  
  return(all_rmp_gene_df)
}

sc_nogene_analysis <- function(hic_data, rmp_data, hic_ct, atac_ct, gencode_granges, signif_th) {
  # filter hic dataset for certain columns, we assume all interactions in this Hi-C dataset are significant
  columns_to_keep <- c('chr1', 'start1', 'end1', 'chr2', 'start2', 'end2', unique(hic_ct))
  hic_data <- hic_data[,columns_to_keep]
  
  # create two sets of granges, one for each set of genomic regions
  hic_granges1 <- GRanges(
    seqnames = hic_data$chr1,
    ranges = IRanges(start = hic_data$start1, end = hic_data$end1),
  )
  
  hic_granges2 <- GRanges(
    seqnames = hic_data$chr2,
    ranges = IRanges(start = hic_data$start2, end = hic_data$end2),
  )
  
  # overlap these granges with promoter granges as defined using gencode earlier
  hic_gene_overlap1 <- findOverlaps(hic_granges1, gencode_granges)
  hic_gene_overlap2 <- findOverlaps(hic_granges2, gencode_granges)
  
  if (length(hic_gene_overlap1) == 0 & length(hic_gene_overlap2) == 0) {
    return('none')
  } else if (length(hic_gene_overlap1) == 0) {
    # create dataframe of results of this overlap
    gene_overlap_df2 <- data.frame(
      Chr1 = hic_data$chr1[queryHits(hic_gene_overlap2)],
      Start1 = hic_data$start1[queryHits(hic_gene_overlap2)],
      End1 = hic_data$end1[queryHits(hic_gene_overlap2)],
      promChr = hic_data$chr2[queryHits(hic_gene_overlap2)],
      promStart = hic_data$start2[queryHits(hic_gene_overlap2)],
      promEnd = hic_data$end2[queryHits(hic_gene_overlap2)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap2)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap2)])
    )
    
    # adding additional information to dataframe
    hic_overlap_df2 <- hic_data[queryHits(hic_gene_overlap2),]
    hic_overlap_df2 <- hic_overlap_df2[,colnames(hic_overlap_df2) %in% hic_ct]
    gene_overlap_df2 <- cbind(gene_overlap_df2, hic_overlap_df2)
    
    # now create granges for other region that interacted with promoter region
    new_hic_granges2 <- GRanges(
      seqnames = gene_overlap_df2$Chr1,
      ranges = IRanges(start = gene_overlap_df2$Start1, end = gene_overlap_df2$End1),
    )
    
    # filter rmp_df based on atac_cell_types
    rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
    
    # create rmp granges
    rmp_granges <- GRanges(
      seqnames = rmp_df_filt$chr,
      ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
    )
    
    hic_rmp_overlap2 <- findOverlaps(new_hic_granges2, rmp_granges)
    
    if (length(hic_rmp_overlap2) == 0) {
      return('none')
    } else {
      # creating dataframe of results of this overlap
      rmp_gene_df2 <- data.frame(
        promChr = gene_overlap_df2$promChr[queryHits(hic_rmp_overlap2)],
        promStart = gene_overlap_df2$promStart[queryHits(hic_rmp_overlap2)],
        promEnd = gene_overlap_df2$promEnd[queryHits(hic_rmp_overlap2)],
        rmpChr = gene_overlap_df2$Chr1[queryHits(hic_rmp_overlap2)],
        rmpStart = gene_overlap_df2$Start1[queryHits(hic_rmp_overlap2)],
        rmpEnd = gene_overlap_df2$End1[queryHits(hic_rmp_overlap2)],
        Gene = gene_overlap_df2$Gene[queryHits(hic_rmp_overlap2)],
        Promoter = gene_overlap_df2$Promoter[queryHits(hic_rmp_overlap2)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap2)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap2)]
      )
      
      # ensuring significance information is still present
      new_overlap_df2 <- gene_overlap_df2[queryHits(hic_rmp_overlap2),]
      new_overlap_df2 <- new_overlap_df2[,colnames(new_overlap_df2) %in% hic_ct]
      rmp_gene_df2 <- cbind(rmp_gene_df2, new_overlap_df2)
      
      all_rmp_gene_df <- rmp_gene_df2 %>% distinct()
    }
  } else if (length(hic_gene_overlap2) == 0) {
    # create dataframe of results of this overlap
    gene_overlap_df1 <- data.frame(
      promChr = hic_data$chr1[queryHits(hic_gene_overlap1)],
      promStart = hic_data$start1[queryHits(hic_gene_overlap1)],
      promEnd = hic_data$end1[queryHits(hic_gene_overlap1)],
      Chr2 = hic_data$chr2[queryHits(hic_gene_overlap1)],
      Start2 = hic_data$start2[queryHits(hic_gene_overlap1)],
      End2 = hic_data$end2[queryHits(hic_gene_overlap1)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap1)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap1)])
    )
    
    # adding additional information from Hi-C dataframe
    hic_overlap_df1 <- hic_data[queryHits(hic_gene_overlap1),]
    hic_overlap_df1 <- hic_overlap_df1[,colnames(hic_overlap_df1) %in% hic_ct]
    gene_overlap_df1 <- cbind(gene_overlap_df1, hic_overlap_df1)
    
    # now create granges for other region that interacted with promoter region
    new_hic_granges1 <- GRanges(
      seqnames = gene_overlap_df1$Chr2,
      ranges = IRanges(start = gene_overlap_df1$Start2, end = gene_overlap_df1$End2),
    )
    
    # filter rmp_df based on atac_cell_types
    rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
    
    # create rmp granges
    rmp_granges <- GRanges(
      seqnames = rmp_df_filt$chr,
      ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
    )
    
    # find overlaps between these granges and RMPs
    hic_rmp_overlap1 <- findOverlaps(new_hic_granges1, rmp_granges)
    
    if (length(hic_rmp_overlap1) == 0) {
      return('none')
    } else {
      # creating dataframe of results of this overlap
      rmp_gene_df1 <- data.frame(
        promChr = gene_overlap_df1$promChr[queryHits(hic_rmp_overlap1)],
        promStart = gene_overlap_df1$promStart[queryHits(hic_rmp_overlap1)],
        promEnd = gene_overlap_df1$promEnd[queryHits(hic_rmp_overlap1)],
        rmpChr = gene_overlap_df1$Chr2[queryHits(hic_rmp_overlap1)],
        rmpStart = gene_overlap_df1$Start2[queryHits(hic_rmp_overlap1)],
        rmpEnd = gene_overlap_df1$End2[queryHits(hic_rmp_overlap1)],
        Gene = gene_overlap_df1$Gene[queryHits(hic_rmp_overlap1)],
        Promoter = gene_overlap_df1$Promoter[queryHits(hic_rmp_overlap1)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap1)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap1)]
      )
      
      # ensuring significance information is still present
      new_overlap_df1 <- gene_overlap_df1[queryHits(hic_rmp_overlap1),]
      new_overlap_df1 <- new_overlap_df1[,colnames(new_overlap_df1) %in% hic_ct]
      rmp_gene_df1 <- cbind(rmp_gene_df1, new_overlap_df1)
      
      all_rmp_gene_df <- rmp_gene_df1 %>% distinct()
    }
  } else {
    # create dataframes of results of this overlap
    gene_overlap_df1 <- data.frame(
      promChr = hic_data$chr1[queryHits(hic_gene_overlap1)],
      promStart = hic_data$start1[queryHits(hic_gene_overlap1)],
      promEnd = hic_data$end1[queryHits(hic_gene_overlap1)],
      Chr2 = hic_data$chr2[queryHits(hic_gene_overlap1)],
      Start2 = hic_data$start2[queryHits(hic_gene_overlap1)],
      End2 = hic_data$end2[queryHits(hic_gene_overlap1)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap1)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap1)])
    )
    
    gene_overlap_df2 <- data.frame(
      Chr1 = hic_data$chr1[queryHits(hic_gene_overlap2)],
      Start1 = hic_data$start1[queryHits(hic_gene_overlap2)],
      End1 = hic_data$end1[queryHits(hic_gene_overlap2)],
      promChr = hic_data$chr2[queryHits(hic_gene_overlap2)],
      promStart = hic_data$start2[queryHits(hic_gene_overlap2)],
      promEnd = hic_data$end2[queryHits(hic_gene_overlap2)],
      Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap2)],
      Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap2)])
    )
    
    # adding additional information from Hi-C dataframe
    hic_overlap_df1 <- hic_data[queryHits(hic_gene_overlap1),]
    hic_overlap_df1 <- hic_overlap_df1[,colnames(hic_overlap_df1) %in% hic_ct]
    gene_overlap_df1 <- cbind(gene_overlap_df1, hic_overlap_df1)
    
    hic_overlap_df2 <- hic_data[queryHits(hic_gene_overlap2),]
    hic_overlap_df2 <- hic_overlap_df2[,colnames(hic_overlap_df2) %in% hic_ct]
    gene_overlap_df2 <- cbind(gene_overlap_df2, hic_overlap_df2)
    
    # now create granges for other region that interacted with promoter region, keeping with the numbering convention of dfs
    new_hic_granges1 <- GRanges(
      seqnames = gene_overlap_df1$Chr2,
      ranges = IRanges(start = gene_overlap_df1$Start2, end = gene_overlap_df1$End2),
    )
    
    new_hic_granges2 <- GRanges(
      seqnames = gene_overlap_df2$Chr1,
      ranges = IRanges(start = gene_overlap_df2$Start1, end = gene_overlap_df2$End1),
    )
    
    # filter rmp_df based on atac_cell_types
    rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
    
    # create rmp granges
    rmp_granges <- GRanges(
      seqnames = rmp_df_filt$chr,
      ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
    )
    
    # find overlaps between these granges and RMPs
    hic_rmp_overlap1 <- findOverlaps(new_hic_granges1, rmp_granges)
    hic_rmp_overlap2 <- findOverlaps(new_hic_granges2, rmp_granges)
    
    # accounting for cases where there is no overlap between RMPs and the filtered Hi-C dataset
    if (length(hic_rmp_overlap1) == 0 & length(hic_rmp_overlap2) == 0) {
      return('none')
    } else if (length(hic_rmp_overlap1) == 0) {
      # creating dataframe of results of this overlap
      rmp_gene_df2 <- data.frame(
        promChr = gene_overlap_df2$promChr[queryHits(hic_rmp_overlap2)],
        promStart = gene_overlap_df2$promStart[queryHits(hic_rmp_overlap2)],
        promEnd = gene_overlap_df2$promEnd[queryHits(hic_rmp_overlap2)],
        rmpChr = gene_overlap_df2$Chr1[queryHits(hic_rmp_overlap2)],
        rmpStart = gene_overlap_df2$Start1[queryHits(hic_rmp_overlap2)],
        rmpEnd = gene_overlap_df2$End1[queryHits(hic_rmp_overlap2)],
        Gene = gene_overlap_df2$Gene[queryHits(hic_rmp_overlap2)],
        Promoter = gene_overlap_df2$Promoter[queryHits(hic_rmp_overlap2)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap2)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap2)]
      )
      
      # ensuring significance information is still present
      new_overlap_df2 <- gene_overlap_df2[queryHits(hic_rmp_overlap2),]
      new_overlap_df2 <- new_overlap_df2[,colnames(new_overlap_df2) %in% hic_ct]
      rmp_gene_df2 <- cbind(rmp_gene_df2, new_overlap_df2)
      
      all_rmp_gene_df <- rmp_gene_df2 %>% distinct()
    } else if (length(hic_rmp_overlap2) == 0) {
      # creating dataframe of results of this overlap
      rmp_gene_df1 <- data.frame(
        promChr = gene_overlap_df1$promChr[queryHits(hic_rmp_overlap1)],
        promStart = gene_overlap_df1$promStart[queryHits(hic_rmp_overlap1)],
        promEnd = gene_overlap_df1$promEnd[queryHits(hic_rmp_overlap1)],
        rmpChr = gene_overlap_df1$Chr2[queryHits(hic_rmp_overlap1)],
        rmpStart = gene_overlap_df1$Start2[queryHits(hic_rmp_overlap1)],
        rmpEnd = gene_overlap_df1$End2[queryHits(hic_rmp_overlap1)],
        Gene = gene_overlap_df1$Gene[queryHits(hic_rmp_overlap1)],
        Promoter = gene_overlap_df1$Promoter[queryHits(hic_rmp_overlap1)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap1)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap1)]
      )
      
      # ensuring significance information is still present
      new_overlap_df1 <- gene_overlap_df1[queryHits(hic_rmp_overlap1),]
      new_overlap_df1 <- new_overlap_df1[,colnames(new_overlap_df1) %in% hic_ct]
      rmp_gene_df1 <- cbind(rmp_gene_df1, new_overlap_df1)
      
      all_rmp_gene_df <- rmp_gene_df1 %>% distinct()
    } else {
      # creating dataframes of results of this overlap
      rmp_gene_df1 <- data.frame(
        promChr = gene_overlap_df1$promChr[queryHits(hic_rmp_overlap1)],
        promStart = gene_overlap_df1$promStart[queryHits(hic_rmp_overlap1)],
        promEnd = gene_overlap_df1$promEnd[queryHits(hic_rmp_overlap1)],
        rmpChr = gene_overlap_df1$Chr2[queryHits(hic_rmp_overlap1)],
        rmpStart = gene_overlap_df1$Start2[queryHits(hic_rmp_overlap1)],
        rmpEnd = gene_overlap_df1$End2[queryHits(hic_rmp_overlap1)],
        Gene = gene_overlap_df1$Gene[queryHits(hic_rmp_overlap1)],
        Promoter = gene_overlap_df1$Promoter[queryHits(hic_rmp_overlap1)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap1)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap1)]
      )
      
      rmp_gene_df2 <- data.frame(
        promChr = gene_overlap_df2$promChr[queryHits(hic_rmp_overlap2)],
        promStart = gene_overlap_df2$promStart[queryHits(hic_rmp_overlap2)],
        promEnd = gene_overlap_df2$promEnd[queryHits(hic_rmp_overlap2)],
        rmpChr = gene_overlap_df2$Chr1[queryHits(hic_rmp_overlap2)],
        rmpStart = gene_overlap_df2$Start1[queryHits(hic_rmp_overlap2)],
        rmpEnd = gene_overlap_df2$End1[queryHits(hic_rmp_overlap2)],
        Gene = gene_overlap_df2$Gene[queryHits(hic_rmp_overlap2)],
        Promoter = gene_overlap_df2$Promoter[queryHits(hic_rmp_overlap2)],
        rmpRegion = GRangesToString(rmp_granges[subjectHits(hic_rmp_overlap2)]),
        cell = rmp_df_filt$cell[subjectHits(hic_rmp_overlap2)]
      )
      
      # ensuring significance information is still present
      new_overlap_df1 <- gene_overlap_df1[queryHits(hic_rmp_overlap1),]
      new_overlap_df1 <- new_overlap_df1[,colnames(new_overlap_df1) %in% hic_ct]
      rmp_gene_df1 <- cbind(rmp_gene_df1, new_overlap_df1)
      
      new_overlap_df2 <- gene_overlap_df2[queryHits(hic_rmp_overlap2),]
      new_overlap_df2 <- new_overlap_df2[,colnames(new_overlap_df2) %in% hic_ct]
      rmp_gene_df2 <- cbind(rmp_gene_df2, new_overlap_df2)
      
      # bind these two dataframes tgt
      all_rmp_gene_df <- rbind(rmp_gene_df1, rmp_gene_df2)
      all_rmp_gene_df <- all_rmp_gene_df %>% distinct()
    }
  }
  
  # filtering out instances for each cell type where Hi-C does not confirm an interaction
  cell_type_combos <- setNames(as.list(hic_ct), atac_ct)
  
  all_results_df <- all_rmp_gene_df
  
  # checking to see if chicago scores or p-values are provided
  if (signif_th >= 1) {
    message("Chicago scores detected as significance values for Hi-C in row ", i, ", filtering for score >= ", signif_th)
    for (cell_type in names(cell_type_combos)) { 
      column_name <- cell_type_combos[[cell_type]]
      all_results_df <- all_results_df %>%
        filter(!(cell == cell_type & get(column_name) < signif_th))
    }
  } else {
    for (cell_type in names(cell_type_combos)) { 
      column_name <- cell_type_combos[[cell_type]]
      all_results_df <- all_results_df %>%
        filter(!(cell == cell_type & get(column_name) > signif_th))
    }
  }
  
  # check to see if all_results_df is now empty
  if (nrow(all_results_df) == 0) {
    return('none')
  } else {
    # finalize this dataframe
    all_results_df$promRegion <- paste(all_results_df$promChr,
                                       all_results_df$promStart, 
                                       all_results_df$promEnd, sep = "-")  
    
    all_results_df$oeRegion <- paste(all_results_df$rmpChr, 
                                     all_results_df$rmpStart, 
                                     all_results_df$rmpEnd, sep = "-")
    
    all_results_df <- all_results_df %>%
      select(-rmpChr, -rmpStart, -rmpEnd, -promChr, -promStart, -promEnd, -Promoter) %>%
      arrange(cell)
    
    colnames(all_results_df)[colnames(all_results_df) == 'Gene'] <- "gene"
    
    all_results_df <- all_results_df[,c('cell', 'gene', 'rmpRegion', 'promRegion', 'oeRegion')] %>% distinct()
    
    return(all_results_df)
  }
}

# when the hic dataset encompasses just one cell type or if it is bulk Hi-C data with gene info
bulk_gene_present_analysis <- function(hic_data, rmp_data, atac_ct) {
  # filter hic dataset for certain columns
  columns_to_keep <- c('baitchr', 'baitstart', 'baitend', 'baitname', 'oechr', 'oestart', 'oeend')
  hic_data <- hic_data[,columns_to_keep]
  
  # create granges for other end (oe) interacting with bait
  hic_granges <- GRanges(
    seqnames = hic_data$oechr,
    ranges = IRanges(start = hic_data$oestart, end = hic_data$oeend),
  )
  
  # filter rmp_df based on atac_cell_types
  rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
  
  # create rmp granges
  rmp_granges <- GRanges(
    seqnames = rmp_df_filt$chr,
    ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
  )
  
  # identify overlaps and create dataframe for results
  overlaps <- findOverlaps(rmp_granges, hic_granges)
  
  # check to ensure one or more overlaps were found
  if (length(overlaps) == 0) {
    return('none')
  } else {
    overlap_df <- data.frame(
      rmp_idx = queryHits(overlaps),
      hic_idx = subjectHits(overlaps),
      rmpChr = as.character(seqnames(rmp_granges)[queryHits(overlaps)]),
      rmpStart = start(rmp_granges)[queryHits(overlaps)],
      rmpEnd = end(rmp_granges)[queryHits(overlaps)],
      cell = rmp_df_filt$cell[queryHits(overlaps)],
      oeChr = as.character(seqnames(hic_granges)[subjectHits(overlaps)]),
      oeStart = start(hic_granges)[subjectHits(overlaps)],
      oeEnd = end(hic_granges)[subjectHits(overlaps)]
    )
    
    # adding additional information from Hi-C dataset
    hic_overlap_df <- hic_data[subjectHits(overlaps),]
    hic_overlap_df <- hic_overlap_df[,!(colnames(hic_overlap_df) %in% c('oechr', 'oestart', 'oeend'))]
    overlap_df <- cbind(overlap_df, hic_overlap_df)
    
    # cleaning up dataframe
    overlap_df_filt <- overlap_df %>%
      mutate(
        rmpRegion = paste(rmpChr, rmpStart, rmpEnd, sep = "-"),
        promRegion = paste(baitchr, baitstart, baitend, sep = "-"),
        oeRegion = paste(oeChr, oeStart, oeEnd, sep = "-")
      ) %>%
      select(-rmpChr, -rmpStart, -rmpEnd, -oeChr, -oeStart, -oeEnd,
             -hic_idx, -rmp_idx, -baitchr, -baitstart, -baitend) %>%
      arrange(cell)
    
    # removing NAs in baitNames
    overlap_df_filt <- overlap_df_filt[!is.na(overlap_df_filt$baitname),] %>% distinct()
    
    # finalizing dataframe
    colnames(overlap_df_filt)[colnames(overlap_df_filt) == 'baitname'] <- "gene"
    
    overlap_df_filt <- overlap_df_filt[,c('cell', 'gene', 'rmpRegion', 'promRegion', 'oeRegion')]
    
    overlap_df_filt <- overlap_df_filt %>%
      separate_rows(gene, sep = ";") %>%
      distinct()
    
    return(overlap_df_filt)
  }
}

sc_gene_present_analysis <- function(hic_data, rmp_data, hic_ct, atac_ct, signif_th) { # when the hic dataset encompasses multiple cell types
  # filter hic dataset for certain columns
  columns_to_keep <- c('baitchr', 'baitstart', 'baitend', 'baitname', 'oechr', 'oestart', 'oeend', unique(hic_ct))
  hic_data <- hic_data[,columns_to_keep]
  
  # create granges for other end (oe) interacting with bait
  hic_granges <- GRanges(
    seqnames = hic_data$oechr,
    ranges = IRanges(start = hic_data$oestart, end = hic_data$oeend),
  )
  
  # filter rmp_df based on atac_cell_types
  rmp_df_filt <- rmp_data[rmp_data$cell %in% atac_ct,]
  
  # create rmp granges
  rmp_granges <- GRanges(
    seqnames = rmp_df_filt$chr,
    ranges = IRanges(start = rmp_df_filt$start, end = rmp_df_filt$end)
  )
  rmp_granges@seqinfo@seqnames <- gsub("^chr", "", rmp_granges@seqinfo@seqnames)

  # identify overlaps and create dataframe for results
  overlaps <- findOverlaps(rmp_granges, hic_granges)
  
  # check to see if there are one or more overlaps
  if (length(overlaps) == 0) {
    return('none')
  } else {
    overlap_df <- data.frame(
      rmp_idx = queryHits(overlaps),
      hic_idx = subjectHits(overlaps),
      rmpChr = as.character(seqnames(rmp_granges)[queryHits(overlaps)]),
      rmpStart = start(rmp_granges)[queryHits(overlaps)],
      rmpEnd = end(rmp_granges)[queryHits(overlaps)],
      cell = rmp_df_filt$cell[queryHits(overlaps)],
      oeChr = as.character(seqnames(hic_granges)[subjectHits(overlaps)]),
      oeStart = start(hic_granges)[subjectHits(overlaps)],
      oeEnd = end(hic_granges)[subjectHits(overlaps)]
    )
    
    # adding additional information from Hi-C dataset
    hic_overlap_df <- hic_data[subjectHits(overlaps),]
    hic_overlap_df <- hic_overlap_df[,!(colnames(hic_overlap_df) %in% c('oechr', 'oestart', 'oeend'))]
    overlap_df <- cbind(overlap_df, hic_overlap_df)
    
    # filtering out instances for each cell type where Hi-C does not confirm an interaction
    cell_type_combos <- setNames(as.list(hic_ct), atac_ct)
    
    overlap_df_filt <- overlap_df
    
    # checking if p-values or chicago scores are provided
    if (signif_th >= 1) {
      message("Chicago scores detected as significance values for Hi-C in row ", i, ", filtering for score >= ", signif_th)
      for (cell_type in names(cell_type_combos)) { 
        column_name <- cell_type_combos[[cell_type]]
        overlap_df_filt <- overlap_df_filt %>%
          filter(!(cell == cell_type & get(column_name) < signif_th))
      }
    } else {
      for (cell_type in names(cell_type_combos)) { 
        column_name <- cell_type_combos[[cell_type]]
        overlap_df_filt <- overlap_df_filt %>%
          filter(!(cell == cell_type & get(column_name) > signif_th))
      }
    }

    # check to see if overlap_df_filt is empty
    if (nrow(overlap_df_filt) == 0) {
      return('none')
    } else {
      # cleaning up dataframe
      overlap_df_filt <- overlap_df_filt %>%
        mutate(
          rmpRegion = paste(rmpChr, rmpStart, rmpEnd, sep = "-"),
          promRegion = paste(baitchr, baitstart, baitend, sep = "-"),
          oeRegion = paste(oeChr, oeStart, oeEnd, sep = "-")
        ) %>%
        select(-rmpChr, -rmpStart, -rmpEnd, -oeChr, -oeStart, -oeEnd,
               -hic_idx, -rmp_idx, -baitchr, -baitstart, -baitend) %>%
        arrange(cell)
      
      # removing NAs in baitNames
      overlap_df_filt <- overlap_df_filt[!is.na(overlap_df_filt$baitname),]
      
      # rename baitName and cell column and make it so there is one gene per row
      colnames(overlap_df_filt)[colnames(overlap_df_filt) == 'baitname'] <- "gene"
      
      overlap_df_filt <- overlap_df_filt[,c('cell', 'gene', 'rmpRegion', 'promRegion', 'oeRegion')]
      
      overlap_df_filt <- overlap_df_filt %>%
        separate_rows(gene, sep = ";") %>%
        distinct()
      
      return(overlap_df_filt)
    }
  }
}

# creating list of results for Hi-C
hic_results_list <- list()

# iterating over number of rows of user_instruct
num_hic <- nrow(user_instruct)

for (i in 1:num_hic) {
  current_row <- user_instruct[i,]
  hic_dir <- current_row$hic_dir
  # read in hic dataset and associated information
  hic_dataset <- read.table(file = hic_dir, sep = '\t', header = TRUE)
  colnames(hic_dataset) = tolower(colnames(hic_dataset))
  genes_present <- as.logical(trimws(current_row$genes_present))
  cell_sorted <- as.logical(trimws(current_row$cell_sorted))
  atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
  
  if (is.na(current_row$hic_cell_types)) {
    hic_cell_types <- NA
  } else {
    hic_cell_types <- unlist(strsplit(current_row$hic_cell_types, split = ","))
  }
  hic_cell_types = tolower(hic_cell_types)
  # run the appropriate function depending on whether genes are present in Hi-C data
  if (genes_present) {
    if (!cell_sorted | is.na(hic_cell_types[1])) { # analysis of bulk is identical to if single cell had just one cell type
      hic_results <- bulk_gene_present_analysis(hic_dataset, rmp_df, atac_cell_types)
      if (is.data.frame(hic_results)) {
        write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                    sep = '\t')
        hic_results_list[[i]] <- hic_results
      } else {
        message("Note: No Hi-C interactions between genes and RMPs found for row with Hi-C directory '", hic_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. One possible reason for this is a mismatch between genome builds of Hi-C, SNP, and GENCODE data.")
      }
    } else { # this means there is one hic dataset encompassing multiple cell types
      signif_threshold <- as.numeric(current_row$signif_threshold)
      hic_results <- sc_gene_present_analysis(hic_dataset, rmp_df, hic_cell_types, atac_cell_types, signif_threshold)
      if (is.data.frame(hic_results)) {
        write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                    sep = '\t')
        hic_results_list[[i]] <- hic_results
      } else {
        message("Note: No Hi-C interactions between genes and RMPs found for row with Hi-C directory '", hic_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. One possible reason for this is a mismatch between genome builds of Hi-C, SNP, and GENCODE data.")
      }
    }
  } else {
    if (!cell_sorted | is.na(hic_cell_types[1])) {
      hic_results <- bulk_nogene_analysis(hic_dataset, rmp_df, atac_cell_types, gene_tss_grg)
      if (is.data.frame(hic_results)) {
        write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                    sep = '\t')
        hic_results_list[[i]] <- hic_results
      } else {
        message("Note: No Hi-C interactions between genes and RMPs found for row with Hi-C directory '", hic_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. One possible reason for this is a mismatch between genome builds of Hi-C, SNP, and GENCODE data.")
      }
    } else { # this means there is one hic dataset encompassing multiple cell types
      signif_threshold <- as.numeric(current_row$signif_threshold)
      hic_results <- sc_nogene_analysis(hic_dataset, rmp_df, hic_cell_types, atac_cell_types, gene_tss_grg, signif_threshold)
      if (is.data.frame(hic_results)) {
        write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                    sep = '\t')
        hic_results_list[[i]] <- hic_results
      } else {
        message("Note: No Hi-C interactions between genes and RMPs found for row with Hi-C directory '", hic_dir, "' ", "and ATAC cell types: '", paste(atac_cell_types, collapse = ','), "' of instructions spreadsheet. One possible reason for this is a mismatch between genome builds of Hi-C, SNP, and GENCODE data.")
      }
    }
  }
}

# check to see if any results were obtained
if (length(hic_results_list) == 0) {
  message('Hi-C analysis complete, no results were found over all Hi-C datasets.')
} else {
  # bind all results together
  all_results <- bind_rows(hic_results_list) %>% distinct()
  rownames(all_results) <- NULL
  
  # save this dataframe
  write.table(all_results, file = paste0(output_dir, "all_hic_results.txt"), row.names = F, quote = F,
              sep = '\t')
}





