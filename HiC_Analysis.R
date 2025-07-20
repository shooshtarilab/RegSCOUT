library(ape)
library(readxl)
library(dplyr)
library(tidyr)
library(readxl)
library(stringr)
library(GenomicRanges)
library(Signac)
library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

# read output directory
output_dir = args[["output_dir"]]

# defining default parameter values
defaults <- list(
  prom_th_up = 2000,
  prom_th_down = 2000
)

# read hic instructions directory
hic_instruct_dir <- args[["hic_instruct_dir"]]

# loading in user instructions
user_instruct <- read_xlsx(hic_instruct_dir)

# load in gencode gene annotation and processing it
prom_th_up = if (!is.null(args[["prom_th_up"]])) as.integer(args[["prom_th_up"]]) else defaults$prom_th_up
prom_th_down = if (!is.null(args[["prom_th_down"]])) as.integer(args[["prom_th_down"]]) else defaults$prom_th_down
gene_annot_dir <- args[["gencode_dir"]]
gene_annot = read.gff(gene_annot_dir, na.strings = c(".", "?"), GFF3 = TRUE)
transcript_type_list = str_split(gene_annot$attributes, "transcript_type=", simplify = TRUE) 
transcript_type_list = str_split(transcript_type_list[,2], ";", simplify = TRUE)[,1]

transcript_coding_index = transcript_type_list == "protein_coding"
gene_annot = gene_annot[transcript_coding_index,]
gene_transcript_data = gene_annot[gene_annot$type == "transcript",]
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

gene_data_temp_pos = gene_transcript_data[gene_transcript_data$strand == "+",]

gene_tss_grg_pos = GRanges(ranges= IRanges(start = gene_data_temp_pos$TSS - prom_th_up,
                                           end = gene_data_temp_pos$TSS + prom_th_down),
                           seqnames = gene_data_temp_pos$seqid)

gene_tss_grg_pos$gene_name = gene_data_temp_pos$gene_name

gene_data_temp_neg = gene_transcript_data[gene_transcript_data$strand == "-",]

gene_tss_grg_neg = GRanges(ranges= IRanges(start = gene_data_temp_neg$TSS - prom_th_down,
                                           end = gene_data_temp_neg$TSS + prom_th_up),
                           seqnames = gene_data_temp_neg$seqid)

gene_tss_grg_neg$gene_name = gene_data_temp_neg$gene_name

gene_tss_grg = c(gene_tss_grg_pos, gene_tss_grg_neg)

# load in RMP information
rmp_df <- read_xlsx(paste0(output_dir, 'risk_regions_ppa.xlsx'))

# create separate chr, start, end columns for RMPs
rmp_df <- rmp_df %>%
  separate(region, into = c("chr", "start", "end"), sep = "-", convert = T)

# expand such that there is only one cell type per row for ease of filtering 
rmp_df <- rmp_df %>%
  separate_rows(cell, sep = ',') %>%
  distinct()

## defining functions for each type of Hi-C analysis (PC-HiC or intact HiC and single cell or bulk)
bulk_intact_analysis <- function(hic_data, rmp_data, atac_ct, gencode_granges) {
  # filter hic dataset for certain columns, we assume all interactions in this Hi-C dataset are significant
  columns_to_keep <- c('Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2')
  hic_data <- hic_data[,columns_to_keep]
  
  # create two sets of granges, one for each set of genomic regions
  hic_granges1 <- GRanges(
    seqnames = hic_data$Chr1,
    ranges = IRanges(start = hic_data$Start1, end = hic_data$End1),
  )
  
  hic_granges2 <- GRanges(
    seqnames = hic_data$Chr2,
    ranges = IRanges(start = hic_data$Start2, end = hic_data$End2),
  )
  
  # overlap these granges with promoter granges as defined using gencode earlier
  hic_gene_overlap1 <- findOverlaps(hic_granges1, gencode_granges)
  hic_gene_overlap2 <- findOverlaps(hic_granges2, gencode_granges)
  
  # create dataframes of results of this overlap
  gene_overlap_df1 <- data.frame(
    promChr = hic_data$Chr1[queryHits(hic_gene_overlap1)],
    promStart = hic_data$Start1[queryHits(hic_gene_overlap1)],
    promEnd = hic_data$End1[queryHits(hic_gene_overlap1)],
    Chr2 = hic_data$Chr2[queryHits(hic_gene_overlap1)],
    Start2 = hic_data$Start2[queryHits(hic_gene_overlap1)],
    End2 = hic_data$End2[queryHits(hic_gene_overlap1)],
    Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap1)],
    Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap1)])
  )
  
  gene_overlap_df2 <- data.frame(
    Chr1 = hic_data$Chr1[queryHits(hic_gene_overlap2)],
    Start1 = hic_data$Start1[queryHits(hic_gene_overlap2)],
    End1 = hic_data$End1[queryHits(hic_gene_overlap2)],
    promChr = hic_data$Chr2[queryHits(hic_gene_overlap2)],
    promStart = hic_data$Start2[queryHits(hic_gene_overlap2)],
    promEnd = hic_data$End2[queryHits(hic_gene_overlap2)],
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
  
  all_rmp_gene_df <- all_rmp_gene_df[,c('cell', 'gene', 'rmpRegion', 'promRegion', 'oeRegion')]
  
  return(all_rmp_gene_df)
}

scintact_analysis <- function(hic_data, rmp_data, hic_ct, atac_ct, gencode_granges, signif_th) {
  # filter hic dataset for certain columns, we assume all interactions in this Hi-C dataset are significant
  columns_to_keep <- c('Chr1', 'Start1', 'End1', 'Chr2', 'Start2', 'End2', unique(hic_ct))
  hic_data <- hic_data[,columns_to_keep]
  
  # create two sets of granges, one for each set of genomic regions
  hic_granges1 <- GRanges(
    seqnames = hic_data$Chr1,
    ranges = IRanges(start = hic_data$Start1, end = hic_data$End1),
  )
  
  hic_granges2 <- GRanges(
    seqnames = hic_data$Chr2,
    ranges = IRanges(start = hic_data$Start2, end = hic_data$End2),
  )
  
  # overlap these granges with promoter granges as defined using gencode earlier
  hic_gene_overlap1 <- findOverlaps(hic_granges1, gencode_granges)
  hic_gene_overlap2 <- findOverlaps(hic_granges2, gencode_granges)
  
  # create dataframes of results of this overlap
  gene_overlap_df1 <- data.frame(
    promChr = hic_data$Chr1[queryHits(hic_gene_overlap1)],
    promStart = hic_data$Start1[queryHits(hic_gene_overlap1)],
    promEnd = hic_data$End1[queryHits(hic_gene_overlap1)],
    Chr2 = hic_data$Chr2[queryHits(hic_gene_overlap1)],
    Start2 = hic_data$Start2[queryHits(hic_gene_overlap1)],
    End2 = hic_data$End2[queryHits(hic_gene_overlap1)],
    Gene = gencode_granges$gene_name[subjectHits(hic_gene_overlap1)],
    Promoter = GRangesToString(gencode_granges[subjectHits(hic_gene_overlap1)])
  )
  
  gene_overlap_df2 <- data.frame(
    Chr1 = hic_data$Chr1[queryHits(hic_gene_overlap2)],
    Start1 = hic_data$Start1[queryHits(hic_gene_overlap2)],
    End1 = hic_data$End1[queryHits(hic_gene_overlap2)],
    promChr = hic_data$Chr2[queryHits(hic_gene_overlap2)],
    promStart = hic_data$Start2[queryHits(hic_gene_overlap2)],
    promEnd = hic_data$End2[queryHits(hic_gene_overlap2)],
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
  
  # filtering out instances for each cell type where Hi-C does not confirm an interaction
  cell_type_combos <- setNames(as.list(hic_ct), atac_ct)
  
  all_results_df <- all_rmp_gene_df
  for (cell_type in names(cell_type_combos)) { 
    column_name <- cell_type_combos[[cell_type]]
    all_results_df <- all_results_df %>%
      filter(!(cell == cell_type & get(column_name) < signif_th))
  }
  
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
  
  all_results_df <- all_results_df[,c('cell', 'gene', 'rmpRegion', 'promRegion', 'oeRegion')]
  
  return(all_results_df)
}

# when the hic dataset encompasses just one cell type or if it is bulk PCHi-C
bulk_promoter_capture_analysis <- function(hic_data, rmp_data, atac_ct) {
  # filter hic dataset for certain columns
  columns_to_keep <- c('baitChr', 'baitStart', 'baitEnd', 'baitName', 'oeChr', 'oeStart', 'oeEnd')
  hic_data <- hic_data[,columns_to_keep]
  
  # create granges for other end (oe) interacting with bait
  hic_granges <- GRanges(
    seqnames = hic_data$oeChr,
    ranges = IRanges(start = hic_data$oeStart, end = hic_data$oeEnd),
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
  hic_overlap_df <- hic_overlap_df[,!(colnames(hic_overlap_df) %in% c('oeChr', 'oeStart', 'oeEnd'))]
  overlap_df <- cbind(overlap_df, hic_overlap_df)
  
  # cleaning up dataframe
  overlap_df_filt <- overlap_df %>%
    mutate(
      rmpRegion = paste(rmpChr, rmpStart, rmpEnd, sep = "-"),
      promRegion = paste(baitChr, baitStart, baitEnd, sep = "-"),
      oeRegion = paste(oeChr, oeStart, oeEnd, sep = "-")
    ) %>%
    select(-rmpChr, -rmpStart, -rmpEnd, -oeChr, -oeStart, -oeEnd,
           -hic_idx, -rmp_idx, -baitChr, -baitStart, -baitEnd) %>%
    arrange(cell)
  
  # removing NAs in baitNames
  overlap_df_filt <- overlap_df_filt[!is.na(overlap_df_filt$baitName),] %>% distinct()
  
  # finalizing dataframe
  colnames(overlap_df_filt)[colnames(overlap_df_filt) == 'baitName'] <- "gene"
  
  overlap_df_filt <- overlap_df_filt %>%
    separate_rows(gene, sep = ";") %>%
    distinct()
  
  return(overlap_df_filt)
}

scpromoter_capture_analysis <- function(hic_data, rmp_data, hic_ct, atac_ct, signif_th) { # when the hic dataset encompasses multiple cell types
  # filter hic dataset for certain columns
  columns_to_keep <- c('baitChr', 'baitStart', 'baitEnd', 'baitName', 'oeChr', 'oeStart', 'oeEnd', unique(hic_ct))
  hic_data <- hic_data[,columns_to_keep]
  
  # create granges for other end (oe) interacting with bait
  hic_granges <- GRanges(
    seqnames = hic_data$oeChr,
    ranges = IRanges(start = hic_data$oeStart, end = hic_data$oeEnd),
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
  hic_overlap_df <- hic_overlap_df[,!(colnames(hic_overlap_df) %in% c('oeChr', 'oeStart', 'oeEnd'))]
  overlap_df <- cbind(overlap_df, hic_overlap_df)
  
  # filtering out instances for each cell type where Hi-C does not confirm an interaction
  cell_type_combos <- setNames(as.list(hic_ct), atac_ct)
  
  overlap_df_filt <- overlap_df
  for (cell_type in names(cell_type_combos)) { 
    column_name <- cell_type_combos[[cell_type]]
    overlap_df_filt <- overlap_df_filt %>%
      filter(!(cell == cell_type & get(column_name) < signif_th))
  }
  
  # cleaning up dataframe
  overlap_df_filt <- overlap_df_filt %>%
    mutate(
      rmpRegion = paste(rmpChr, rmpStart, rmpEnd, sep = "-"),
      promRegion = paste(baitChr, baitStart, baitEnd, sep = "-"),
      oeRegion = paste(oeChr, oeStart, oeEnd, sep = "-")
    ) %>%
    select(-rmpChr, -rmpStart, -rmpEnd, -oeChr, -oeStart, -oeEnd,
           -hic_idx, -rmp_idx, -baitChr, -baitStart, -baitEnd, -tCD8, -tCD4, 
           -tB) %>%
    arrange(cell)
  
  # removing NAs in baitNames
  overlap_df_filt <- overlap_df_filt[!is.na(overlap_df_filt$baitName),]
  
  # rename baitName and cell column and make it so there is one gene per row
  colnames(overlap_df_filt)[colnames(overlap_df_filt) == 'baitName'] <- "gene"
  
  overlap_df_filt <- overlap_df_filt %>%
    separate_rows(gene, sep = ";") %>%
    distinct()
  
  return(overlap_df_filt)
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
  hic_type <- current_row$hic_type
  bulk <- as.logical(str_to_upper(current_row$bulk))
  atac_cell_types <- unlist(strsplit(current_row$atac_cell_types, split = ","))
  
  if (is.na(current_row$hic_cell_types)) {
    hic_cell_types <- NA
  } else {
    hic_cell_types <- unlist(strsplit(current_row$hic_cell_types, split = ","))
  }
  
  # run the appropriate function depending on Hi-C dataset
  if (hic_type == "promoter capture") {
    if (bulk | is.na(hic_cell_types[1])) { # analysis of bulk is identical to if scPCHi-C had just one cell type
      hic_results <- bulk_promoter_capture_analysis(hic_dataset, rmp_df, atac_cell_types)
      write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
      hic_results_list[[i]] <- hic_results
    } else { # this means there is one hic dataset encompassing multiple cell types
      signif_threshold <- current_row$signif_threshold
      hic_results <- scpromoter_capture_analysis(hic_dataset, rmp_df, hic_cell_types, atac_cell_types, signif_threshold)
      write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
      hic_results_list[[i]] <- hic_results
    }
  } else if (hic_type == "intact") {
    if (bulk | is.na(hic_cell_types[1])) {
      hic_results <- bulk_intact_analysis(hic_dataset, rmp_df, atac_cell_types, gene_tss_grg)
      write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
      hic_results_list[[i]] <- hic_results
    } else { # this means there is one hic dataset encompassing multiple cell types
      signif_threshold <- current_row$signif_threshold
      hic_results <- scintact_analysis(hic_dataset, rmp_df, hic_cell_types, atac_cell_types, gene_tss_grg, signif_threshold)
      write.table(hic_results, file = paste0(output_dir, "hic_analysis", i, "_results.txt"), row.names = F, quote = F,
                  sep = '\t')
      hic_results_list[[i]] <- hic_results
    }
  }
}

# bind all results together
all_results <- bind_rows(hic_results_list)
rownames(all_results) <- NULL

# save this dataframe
write.table(all_results, file = paste0(output_dir, "all_hic_results.txt"), row.names = F, quote = F,
            sep = '\t')


print('Hi-C analysis complete!')
