suppressPackageStartupMessages(library(atSNP))
suppressPackageStartupMessages(library(TFBSTools))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(R.utils))

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

#Defining default parameter values
defaults <- list(
  jaspar_mtx_dir = "none"
)

# read output directory
output_dir = args[["output_dir"]]

# read in and prepare JASPAR matrix file
jaspar_mtx_dir <- if (nzchar(args[["jaspar_mtx_dir"]])) {
  args[["jaspar_mtx_dir"]]
} else {
  message("Using default jaspar_mtx_dir value: ", defaults$jaspar_mtx_dir)
  defaults$jaspar_mtx_dir
}

if (jaspar_mtx_dir == 'none') {
  # load required libraries for obtaining JASPAR matrix
  suppressPackageStartupMessages(library(JASPAR2024))
  
  jaspar_sql = JASPAR2024()
  
  jaspar_pwm <- getMatrixSet(
    x = jaspar_sql@db,
    opts = list(species = 9606, all_versions = FALSE)
  )
} else {
  jaspar_pwm <- readJASPARMatrix(jaspar_mtx_dir, matrixClass = c("PFM", "PWM", "PWMProb"))
}

jaspar_pwm_list = jaspar_pwm@listData

jaspar_pwm_atsnp = list()

for (i in c(1:length(jaspar_pwm_list))){
  temp_tf = jaspar_pwm_list[[i]]
  temp_tf_name = temp_tf@name
  temp_tf_matrix = temp_tf@profileMatrix
  new_matrix = apply(temp_tf_matrix,2,function(x)(x/sum(x)))
  jaspar_pwm_atsnp[[temp_tf_name]] = t(new_matrix)
}

# reading in GWAS data
ci_gwas_dir = args[["ci_gwas_dir"]]
ci_gwas_data = read.table(file=ci_gwas_dir, sep="", header=TRUE)

# Preparing CI SNPs for motif analysis
genome_build = args[["genome_build"]]

if (genome_build == "hg19"){
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg19))
  eff_snp = ci_gwas_data
  
  snp_table = as.data.frame(matrix(0, nrow = nrow(eff_snp), ncol = 5))
  colnames(snp_table) = c("chr","snp","snpid","a1","a2")
  snp_table$chr = eff_snp$chr
  snp_table$snp = eff_snp$pos
  snp_table$snpid = eff_snp$id
  snp_table$a1 = eff_snp$a1
  snp_table$a2 = eff_snp$a2
  
  snp_table_dir = paste0(getwd(), "/", "snp_table.txt")
  write.table(snp_table, file = snp_table_dir, col.names = TRUE, 
              row.names = FALSE, quote = FALSE)
  
  reg1_snp_data = LoadSNPData(filename = snp_table_dir,
                              genome.lib ="BSgenome.Hsapiens.UCSC.hg19"
                              , half.window.size = 30, default.par = TRUE
                              , mutation = FALSE)
  invisible(file.remove(snp_table_dir))
}else if(genome_build == "hg38"){
  suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
  
  eff_snp = ci_gwas_data
  
  snp_table = as.data.frame(matrix(0, nrow = nrow(eff_snp), ncol = 5))
  colnames(snp_table) = c("chr","snp","snpid","a1","a2")
  snp_table$chr = eff_snp$chr
  snp_table$snp = eff_snp$pos
  snp_table$snpid = eff_snp$id
  snp_table$a1 = eff_snp$a1
  snp_table$a2 = eff_snp$a2
  
  snp_table_dir = paste0(getwd(), "/", "snp_table.txt")
  write.table(snp_table, file = snp_table_dir, col.names = TRUE, 
              row.names = FALSE, quote = FALSE)
  
  reg1_snp_data = LoadSNPData(filename = snp_table_dir,
                              genome.lib ="BSgenome.Hsapiens.UCSC.hg38"
                              , half.window.size = 30, default.par = TRUE
                              , mutation = FALSE)
  invisible(file.remove(snp_table_dir))
}

# Calculating Effect SNPs
results = suppressMessages(ComputeMotifScore(jaspar_pwm_atsnp, reg1_snp_data, ncores = 2))

# Initialize an empty list to store the results
all_results <- list()

# Loop ComputePValues() function 10 times and obtain all results
for (i in 1:10) {
  results_pval_i <- suppressMessages((ComputePValues(jaspar_pwm_atsnp, reg1_snp_data, results$motif.scores, ncores = 2, testing.mc=T)))
  all_results[[i]] <- results_pval_i
  message(paste0("Running ComputePValues step ",i)) # progress bar otherwise this can take too long and can look like the pipeline is stuck
} 

# Combine all data frames into one big data frame
results_pval <- do.call(rbind, all_results)

# save this dataframe as RDS
saveRDS(results_pval, file = paste0(output_dir,'atSNP_10runs_results.RDS'))

# Correction for multiple testing
results_pval_val = results_pval$pval_diff
results_pval_val_cor = p.adjust(results_pval_val, method = "fdr")
results_pval$results_pval_val_cor <- results_pval_val_cor

data_with_counts <- results_pval %>%
  group_by(snpid, motif) %>%
  mutate(
    Count = n(),  # total count of each SNP-TF pair
    Count_FDR_lt_0.05 = sum(results_pval_val_cor < 0.05)  # count of FDR values < 0.05
  )

# identifying SNP-TF pairs found to be significant in atleast two iterations of ComputePValue
results_pval_sig = data_with_counts[data_with_counts$Count_FDR_lt_0.05 >= 2,]
results_pval_sig = results_pval_sig[results_pval_sig$results_pval_val_cor < 0.05,]

# creating final dataframe
reg1_effect_snps = sort(results_pval_sig$snpid)

final_effect_snp_frame = data.frame(matrix(data = NA, 
                                           nrow = length(results_pval_sig$snpid),
                                           ncol = 9))
colnames(final_effect_snp_frame) = c("SNP","CHR","Pos","Locus", "TF", "Raw_P_value", "FDR_corr_P_value",
                                     "log_like_ratio", "ppa")

for (i in c(1:length(results_pval_sig$snpid))){
  temp_id = results_pval_sig$snpid[i]
  ci_effect_index = which(ci_gwas_data$id == temp_id)
  temp_pos = ci_gwas_data$pos[ci_effect_index]
  temp_chr = ci_gwas_data$chr[ci_effect_index]
  temp_locus = ci_gwas_data$chunk[ci_effect_index]
  temp_TF <- gsub("^MA[0-9]+\\.[0-9]+\\.", "", results_pval_sig$motif[i])
  temp_pval = results_pval_sig$pval_diff[i]
  temp_pval_corr = results_pval_sig$results_pval_val_cor[i]
  temp_log_ratio = results_pval_sig$log_lik_ratio[i]
  
  final_effect_snp_frame[i,]$ppa = ci_gwas_data$PPA[ci_effect_index]
  final_effect_snp_frame[i,]$SNP = temp_id
  final_effect_snp_frame[i,]$CHR = temp_chr
  final_effect_snp_frame[i,]$Pos = temp_pos
  final_effect_snp_frame[i,]$Locus = temp_locus
  final_effect_snp_frame[i,]$TF = temp_TF
  final_effect_snp_frame[i,]$FDR_corr_P_value = temp_pval_corr
  final_effect_snp_frame[i,]$Raw_P_value = temp_pval
  final_effect_snp_frame[i,]$log_like_ratio = temp_log_ratio
}

Ci_effect_SNPs <- final_effect_snp_frame %>%
  distinct() %>%
  group_by(SNP, CHR, Pos, Locus, TF, ppa) %>%
  summarize(
    FDR_corr_P_value = min(FDR_corr_P_value, na.rm = TRUE),  
    Raw_P_value = min(Raw_P_value, na.rm = TRUE),  
    log_like_ratio = log_like_ratio[which.max(abs(log_like_ratio))],
    .groups = "drop"
  )

Ci_effect_SNPs$log_like_ratio <- -(Ci_effect_SNPs$log_like_ratio) 

final_ci_effect_dir = paste0(output_dir,"Ci_effect_SNPs.txt")
write.table(file = final_ci_effect_dir, Ci_effect_SNPs, col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")


message("Effect-SNP identification complete!")
