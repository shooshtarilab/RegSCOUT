library(GenomicRanges)
#library(liftOver)
#library(rtracklayer)
#library(Repitools)
#library(LDlinkR)
library(dplyr)
#library(ggplot2)
#library(qqman)
#library(TwoSampleMR)
library(genetics.binaRies)
#library(bigsnpr)
library(R.utils)

ld_matrix_local_mod <- function(variants, bfile, plink_bin, fn_file, with_alleles=FALSE)
{
  # Make textfile
  shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
  #fn <- "/Users/naderhosseininaghavi/rprojects_whole/ms_pbmc/gwas_data/bmi_chip/ref_data/temp"
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


#ld_mat = ld_matrix_local(
#  snp_list,
#  plink_bin = genetics.binaRies::get_plink_binary(),
#  bfile = "/Users/naderhosseininaghavi/rprojects_whole/ms_pbmc/gwas_data/bmi_chip/ref_data/g1000_eur"
#)


args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Setting the output directory
#output_dir = args[["output_dir"]]
output_dir = "/Users/naderhosseininaghavi/rprojects_whole/Final_new_pip/results/MS/10x/"

#Getting the Plink files for reference SNPs
#SNP_ref = args[["SNP_ref"]]
SNP_ref = "/Users/naderhosseininaghavi/rprojects_whole/Final_new_pip/results/MS/10x/1kg.v3/"

#Getting the population 
#snp_population = args[["Population"]]
snp_population = "EUR"

snp_ref_files = paste0(SNP_ref,snp_population)

#Loading the summary statistics file
#gwas_data_dir = args[["sum_stats"]]
gwas_data_dir = "/Users/naderhosseininaghavi/rprojects_whole/Final_new_pip/ms_data/new_discovery"
gwas_data = read.table(gwas_data_dir, sep = "\t", header = TRUE)
gwas_data = gwas_data[!is.na(gwas_data$P),]
gwas_data = gwas_data[gwas_data$P < 1e-02,]

gwas_data$CHR = gsub("chr", "", tolower(gwas_data$CHR))


#Loading the lead SNP data
#loci_head_dir = args[["lead_snps"]]
loci_head_dir = "/Users/naderhosseininaghavi/rprojects_whole/ms_pbmc/gwas_data/bmi_chip/lead_SNPs_filt.txt"
loci_head = read.table(loci_head_dir, sep = "\t", header = TRUE)


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
  #Getting the allele frequeny file from reference panel SNP data using Plink
  freq_file_out = paste0(output_dir,"Plink2")
  #Getting the directory of Plink2 software
  #plink2_bin = args[["plink2_bin"]]
  plink2_bin = "/Users/naderhosseininaghavi/rprojects_whole/Final_new_pip/code/plink2"
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
new_loci_dir = paste0(output_dir,"new_loci.txt")
write.table(new_loci_head, file = new_loci_dir, col.names = TRUE,
            row.names = FALSE, quote = FALSE)

#Defining loci based on lead SNPs and assigning SNPs to the loci
for (i in c(1:length(new_loci_head$SNP))){
  #Finding the interval and chr of a locus based on its lead SNP
  index = which(new_gwas_data$SNP == new_loci_head$SNP[i])
  interval = c((new_gwas_data$POS[index]- 1e6):(new_gwas_data$POS[index] + 1e+6))
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
  found_index = abs(ld_signal) > 0.25
  snp_list_loc = snp_list[found_index]
  
  #Assigning locus id to the SNPs
  new_gwas_data[(new_gwas_data$SNP %in% snp_list_loc),'SEGNUM'] = i
  print(paste0("step ",i," completed"))
  #Sys.sleep(60) 
}

print("SNPs assigned to loci!")

#Filtering the SNPs outside of the loci
filt_gwas_data = new_gwas_data[new_gwas_data$SEGNUM != 'None',]

#Creating the final GWAS dataframe 
final_gwas_data = as.data.frame(matrix(nrow = nrow(filt_gwas_data), ncol = 9))
colnames(final_gwas_data) = c("SNPID", "CHR", "POS", "Z", "F", "N", "SEGNUMBER",
                              "A1","A2")

#Setting the sample number
#sample_number = args[["sample_num"]]
sample_number = 115803

final_gwas_data$SNPID = filt_gwas_data$SNP
final_gwas_data$CHR = paste0("chr", filt_gwas_data$CHR)
final_gwas_data$POS = filt_gwas_data$POS
final_gwas_data$Z = filt_gwas_data$Z
final_gwas_data$F = filt_gwas_data$MAF
final_gwas_data$N = rep(sample_number, nrow(final_gwas_data)) 
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

