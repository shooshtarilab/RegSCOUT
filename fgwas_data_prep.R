library(GenomicRanges)
library(dplyr)
library(genetics.binaRies)
library(R.utils)

ld_matrix_local_mod <- function(variants, bfile, plink_bin, fn_file, with_alleles=FALSE)
{
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  fn = fn_file
  write.table(data.frame(variants), file=fn, row.names=F, col.names=F, quote=F)
  
  
  fun1 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --extract ", shQuote(fn, type=shell), 
    " --make-just-bim ", 
    " --keep-allele-order ",
    " --out ", shQuote(fn, type=shell)
  )
  system(fun1)
  
  bim <- read.table(paste0(fn, ".bim"), stringsAsFactors=FALSE)
  
  fun2 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --extract ", shQuote(fn, type=shell), 
    " --r square ", 
    " --keep-allele-order ",
    " --out ", shQuote(fn, type=shell)
  )
  system(fun2)
  
  fun3 <- paste0(
    shQuote(plink_bin, type=shell),
    " --bfile ", shQuote(bfile, type=shell),
    " --extract ", shQuote(fn, type=shell), 
    " --freqx ", 
    " --keep-allele-order ",
    " --out ", shQuote(fn, type=shell)
  )
  system(fun3)
  
  
  res <- read.table(paste0(fn, ".ld"), header=FALSE) %>% as.matrix
  if(with_alleles)
  {
    rownames(res) <- colnames(res) <- paste(bim$V2, bim$V5, bim$V6, sep="_")
  } else {
    rownames(res) <- colnames(res) <- bim$V2
  }
  return(res)
}

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Defining default parameter values
defaults <- list(
  sample_num = "present",
  locus_region = 1000000,
  ld_th = 0.25
)

#Setting the output directory
output_dir = args[["output_dir"]]

#Getting the Plink files for reference SNPs
snp_ref = args[["snp_ref_dir"]]

#Getting the population 
snp_population = args[["population"]]

snp_ref_files = paste0(snp_ref,snp_population)

#Loading the summary statistics file
gwas_data_dir = args[["sum_stats_dir"]]
gwas_data = read.table(gwas_data_dir, sep = "\t", header = TRUE)
gwas_data = gwas_data[!is.na(gwas_data$P),]
gwas_data = gwas_data[gwas_data$P < 1e-02,]

gwas_data$CHR = gsub("chr", "", tolower(gwas_data$CHR))

#Loading the lead SNP data
loci_head_dir = args[["lead_snps_dir"]]
loci_head = read.table(loci_head_dir, sep = "\t", header = TRUE)

# defining a locus number column to be filled in
gwas_data['SEGNUM'] = 'None'

#Caclulating the Z score value for SNPs based on the P values
gwas_qval = qnorm(gwas_data$P, mean = 0, sd = 1, lower.tail = FALSE)
gwas_data['Z'] = gwas_qval

snp_list = gwas_data$SNP

if (sum(c("A1","A2","MAF") %in% colnames(gwas_data)) == 3){
  new_gwas_data = gwas_data
  new_gwas_data[,"A1"] = substr(gwas_data$A1,1,1)
  new_gwas_data[,"A2"] = substr(gwas_data$A2,1,1)
  new_gwas_data[,"MAF"] = gwas_data$MAF
}else{
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  #Getting the allele frequency file from reference panel SNP data using Plink
  freq_file_out = paste0(output_dir,"Plink2")
  #Getting the directory of Plink2 software
  plink2_bin = args[["plink2_bin"]]
  fun <- paste0(
    shQuote(plink2_bin, type=shell),
    " --bfile ", shQuote(snp_ref_files, type=shell),
    " --freq ", 
    " --out ", shQuote(freq_file_out, type=shell)
  )
  system(fun)
  
  #Loading the Plink afreq reference file for SNPs
  snp_freq_file = paste0(freq_file_out,".afreq")

  snp_freq = read.table(snp_freq_file, sep = "\t", header = FALSE)
  unique_snp_freq = snp_freq %>% group_by(V2) %>% filter(row_number() == 1)
  
  #Matching the given SNPs to the reference SNPs
  matched_index = intersect(gwas_data$SNP, unique_snp_freq$V2)
  new_gwas_data = gwas_data[gwas_data$SNP %in% matched_index,]
  unique_snp_freq = unique_snp_freq[unique_snp_freq$V2 %in% matched_index,]
  new_gwas_data = new_gwas_data[match(unique_snp_freq$V2, new_gwas_data$SNP),]
  new_gwas_data[,"A1"] = substr(unique_snp_freq$V3,1,1)
  new_gwas_data[,"A2"] = substr(unique_snp_freq$V4,1,1)
  
  new_gwas_data["MAF"] = unique_snp_freq$V6
}
print("SNP information loaded!")

#Matching lead SNPs to the GWAS SNPs
new_loci_head = loci_head[loci_head$SNP %in% new_gwas_data$SNP,]
#Saving the matched lead SNPs
new_loci_dir = paste0(output_dir,"new_lead_snps.txt")
write.table(new_loci_head, file = new_loci_dir, col.names = TRUE,
            row.names = FALSE, quote = FALSE)

#Defining loci based on lead SNPs and assigning SNPs to the loci
new_gwas_data$lead_snp <- 'None'
new_gwas_data$locus_chr <- NA
new_gwas_data$locus_start <- NA
new_gwas_data$locus_end <- NA

locus_region <- if (nzchar(args[["locus_region"]])) {
  as.integer(args[["locus_region"]])
} else {
  message("Using default locus_region value: ", defaults$locus_region)
  defaults$locus_region
}

ld_th <- if (nzchar(args[["ld_th"]])) {
  as.numeric(args[["ld_th"]])
} else {
  message("Using default ld_th value: ", defaults$ld_th)
  defaults$ld_th
}

seg_index <- 0

for (i in c(1:length(new_loci_head$SNP))){
  #Finding the interval and chr of a locus based on its lead SNP
  index = which(new_gwas_data$SNP == new_loci_head$SNP[i])
  interval = c((new_gwas_data$POS[index]- locus_region):(new_gwas_data$POS[index] + locus_region))
  index_chr = new_gwas_data$CHR[index]
  
  #Finding the SNPs falling into the locus
  found_index = which((new_gwas_data$CHR == index_chr) & (new_gwas_data$POS %in% interval))
  snp_list = new_gwas_data$SNP[found_index]
  
  #Setting population reference files path and prefix
  fn_temp = paste0(output_dir,"temp")
  #Calculating the ld matrix in the locus
  ld_mat = ld_matrix_local_mod(
    snp_list,
    plink_bin = genetics.binaRies::get_plink_binary(),
    bfile = snp_ref_files,
    fn_file = fn_temp
  )
  
  #Only including the SNPs found in the ld data
  ld_found_index = snp_list %in% colnames(ld_mat)
  snp_list = snp_list[ld_found_index]
  new_gwas_loc = new_gwas_data[new_gwas_data$SNP %in% snp_list,]
  if (new_loci_head$SNP[i] %in% colnames(ld_mat)){
    new_lead_snp = new_loci_head$SNP[i]
  }else{
    new_lead_snp = new_gwas_loc$SNP[new_gwas_loc$Z == max(new_gwas_loc$Z)] 
  }
  
  #Filtering the SNPs with ld value of less than 0.25 with the lead SNP
  ld_signal = ld_mat[,new_lead_snp]
  found_index = abs(ld_signal) > ld_th
  snp_list_loc = snp_list[found_index]
  
  #Assigning locus id to the SNPs
  if (length(snp_list_loc) > 0) {
    new_gwas_data[(new_gwas_data$SNP %in% snp_list_loc),'SEGNUM'] = seg_index
    seg_index = seg_index + 1 # ensuring segnumber values always start at zero
    
    #Providing lead SNP and locus region information
    new_gwas_data[(new_gwas_data$SNP %in% snp_list_loc),'lead_snp'] = new_lead_snp
    new_gwas_data[(new_gwas_data$SNP %in% snp_list_loc),'locus_chr'] = index_chr
    new_gwas_data[(new_gwas_data$SNP %in% snp_list_loc),'locus_start'] = new_gwas_data$POS[index] - locus_region
    new_gwas_data[(new_gwas_data$SNP %in% snp_list_loc),'locus_end'] = new_gwas_data$POS[index] + locus_region
  }
  
  print(paste0("step ",i," completed"))
}

print("SNPs assigned to loci!")

#Filtering the SNPs outside of the loci
filt_gwas_data = new_gwas_data[new_gwas_data$SEGNUM != 'None',]

#Creating a file with just locus number and region information
loci_df = filt_gwas_data[, c('SEGNUM', 'lead_snp', 'locus_chr', 'locus_start', 'locus_end')] %>% distinct()
loci_df$SEGNUM <- as.integer(loci_df$SEGNUM)
loci_df <- loci_df[order(loci_df$SEGNUM),]
colnames(loci_df)[colnames(loci_df) == "SEGNUM"] <- "chunk"
loci_df$locus_start[loci_df$locus_start < 0] <- 0
loci_df$locus_chr = paste0("chr", loci_df$locus_chr)

loci_info_dir = paste0(output_dir,"loci_info.txt")
write.table(loci_df, file = loci_info_dir, col.names = TRUE, sep = '\t',
            row.names = FALSE, quote = FALSE)

#Creating the final GWAS dataframe 
final_gwas_data = as.data.frame(matrix(nrow = nrow(filt_gwas_data), ncol = 9))
colnames(final_gwas_data) = c("SNPID", "CHR", "POS", "Z", "F", "N", "SEGNUMBER",
                              "A1","A2")

#Getting value of sample_num flag
sample_num <- if (nzchar(args[["sample_num"]])) {
  args[["sample_num"]]
} else {
  message("Using default sample_num value: ", defaults$sample_num)
  defaults$sample_num
}

final_gwas_data$SNPID = filt_gwas_data$SNP
final_gwas_data$CHR = paste0("chr", filt_gwas_data$CHR)
final_gwas_data$POS = filt_gwas_data$POS
final_gwas_data$Z = filt_gwas_data$Z
final_gwas_data$F = filt_gwas_data$MAF
if (sample_num == "present") {
  final_gwas_data$N = filt_gwas_data$N
} else {
  final_gwas_data$N = rep(as.integer(sample_num), nrow(final_gwas_data)) # setting a constant sample number for each SNP, if N column not present
}
final_gwas_data$SEGNUMBER = filt_gwas_data$SEGNUM
final_gwas_data$A1 = filt_gwas_data$A1
final_gwas_data$A2 = filt_gwas_data$A2

#Setting the output file
final_gwas_dir = paste0(output_dir,"final_gwas_data.txt")
final_gwas_data = final_gwas_data[order(as.numeric(final_gwas_data$SEGNUMBER)),]

#Ordering the data based on locus number
for (j in c(1:i)){
  found_index = final_gwas_data$SEGNUMBER == j
  found_data = final_gwas_data[found_index,]
  final_gwas_data[found_index,] = found_data[order(found_data$POS),]
  print(j)
}

#Saving the final table
final_gwas_data <- final_gwas_data[!is.infinite(final_gwas_data$Z),]
write.table(final_gwas_data, file = final_gwas_dir, sep = "\t", col.names = TRUE,
            row.names = FALSE, quote = FALSE)

print("SNPs are ready for fine mapping!")

