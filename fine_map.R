suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#defining defaults for parameters
defaults <- list(
  ci_th = 0.95,
  ci_ppa_th = 0.01
)

output_dir = args[["output_dir"]]

fgwas_src = args[["fgwas_dir"]]

if (!nzchar(fgwas_src)){
    fgwas_src = Sys.which("fgwas")
}

fgwas_file = paste0(output_dir,"final_gwas_data.txt")
loci_info_file = paste0(output_dir, "loci_info_temp.txt")
ci_suff = "CI"
ci_files = paste0(output_dir,ci_suff)

shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")

fun <- paste0(
    shQuote(fgwas_src, type=shell),
    " -i ", shQuote(fgwas_file, type=shell),
    " -o ", shQuote(ci_files, type=shell),
    " -fine ",
    " -print"
)
system(fun)

ci_file_bfs = paste0(ci_files,".bfs.gz")

fun <- paste0(
    shQuote("gzip", type=shell),
    " -d ",
    " -f ", # force overwrite a problem?
    shQuote(ci_file_bfs, type=shell)
)

system(fun)

ci_files_bfs = paste0(ci_files,".bfs")
ci_data = read.table(ci_files_bfs, sep = " ", header = TRUE)
fgwas_data = read.table(fgwas_file, sep = "\t", header = TRUE)
loci_info = read.table(loci_info_file, sep = "\t", header = TRUE)

# adding loci information to fgwas data
fgwas_data <- fgwas_data %>% 
  left_join(loci_info, by = c("SEGNUMBER" = "chunk"))

fgwas_data <- fgwas_data %>% 
  mutate(SEGNUMBER = dense_rank(SEGNUMBER))

# outputting new loci_info text file, deleting previous file
loci_df = fgwas_data[, c('SEGNUMBER', 'lead_snp', 'locus_chr', 'locus_start', 'locus_end')] %>% distinct()
loci_df$SEGNUMBER <- as.integer(loci_df$SEGNUMBER)
loci_df <- loci_df[order(loci_df$SEGNUMBER),]
colnames(loci_df)[colnames(loci_df) == "SEGNUMBER"] <- "chunk"

loci_info_dir = paste0(output_dir,"loci_info.txt")
write.table(loci_df, file = loci_info_dir, col.names = TRUE, sep = '\t',
            row.names = FALSE, quote = FALSE)

unlink(loci_info_file)

# outputting new final_gwas_data text file, overwriting previous
columns_to_remove <- c('lead_snp', 'locus_chr', 'locus_start', 'locus_end')
fgwas_data_out <- fgwas_data[,!(colnames(fgwas_data) %in% columns_to_remove)]

final_gwas_dir = paste0(output_dir,"final_gwas_data.txt")
write.table(fgwas_data_out, file = final_gwas_dir, sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

# ensuring locus/chunk numbers remain consistent
ci_data$locus <- NA
ci_data$locus <- fgwas_data$SEGNUMBER[match(ci_data$id, fgwas_data$SNPID)] 
ci_data$chunk <- ci_data$locus
ci_data$locus <- NULL

# adding reference and alternate allele information to fgwas results
ci_data$a1 <- fgwas_data$A1[match(ci_data$id, fgwas_data$SNPID)] 
ci_data$a2 <- fgwas_data$A2[match(ci_data$id, fgwas_data$SNPID)]

## Removing SNPs with more than one base as either the alternate or reference allele
ci_data = ci_data[which(nchar(as.character(ci_data$a1)) == 1 & nchar(as.character(ci_data$a2)) == 1),]

# obtaining smallest number of SNPs whose PPAs sum up to ci_th
ci_th <- if (is.null(args[["ci_th"]])) {
  message("Using default ci_th value: ", defaults$ci_th)
  defaults$ci_th
} else if (!nzchar(args[["ci_th"]])) {
  message("Using default ci_th value: ", defaults$ci_th)
  defaults$ci_th
} else {
  as.numeric(args[["ci_th"]])
}

ci_final_list = list()
region_list = unique(ci_data$chunk)
for (i in region_list){
  region_data = ci_data[ci_data$chunk == i,]
  region_data = region_data[order(region_data$PPA, decreasing = TRUE),]
  ci_index = c()
  ci_sum = 0
  count = 1
  while (ci_sum < ci_th){
    ci_index = append(ci_index, count)
    ci_sum = ci_sum + region_data[count,]$PPA
    count = count + 1
    if (count > nrow(region_data)){
      break
    }
  }
  ci_final_list = append(ci_final_list, region_data[ci_index,]$id)
}
ci_final_list = unlist(ci_final_list)
ci_gwas_data = ci_data[ci_data$id %in% ci_final_list,]

# removing all SNPs on sex chromosomes
chr_list = paste0("chr", c(1:22))
ci_gwas_data = ci_gwas_data[ci_gwas_data$chr %in% chr_list,]

# removing all SNPs with PPA <= CI PPA threshold
ci_ppa_th = if (is.null(args[["ci_ppa_th"]])) {
  message("Using default ci_ppa_th value: ", defaults$ci_th)
  defaults$ci_th
} else if (!nzchar(args[["ci_ppa_th"]])) {
  message("Using default ci_ppa_th value: ", defaults$ci_ppa_th)
  defaults$ci_ppa_th
} else {
  as.numeric(args[["ci_ppa_th"]])
} 

ci_gwas_data = ci_gwas_data[ci_gwas_data$PPA > ci_ppa_th,]

ci_dir = paste0(output_dir,"gwas_CI.txt")
write.table(ci_gwas_data, file = ci_dir, col.names = TRUE, sep="\t",
            row.names = FALSE, quote = FALSE)

files_to_delete <- file.path(output_dir, c("CI.bfs", "CI.llk", "CI.params", "CI.ridgeparams", "CI.segbfs.gz"))

# delete the files
invisible(suppressWarnings(file.remove(files_to_delete)))
