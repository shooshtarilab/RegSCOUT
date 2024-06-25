library(atSNP)
library(JASPAR2024)
library(TFBSTools)
library(R.utils)
library(GenomicRanges)
library(dplyr)
library(RSQLite)
#library(SNPlocs.Hsapiens.dbSNP144.GRCh37)
#library(SNPlocs.Hsapiens.dbSNP144.GRCh38)


args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

jaspar_sql = JASPAR2024()

jaspar_pwm <- getMatrixSet(
  x = jaspar_sql@db,
  opts = list(species = 9606, all_versions = FALSE)
)

print("JASPAR PWM matrix loaded!")


jaspar_pwm_list = jaspar_pwm@listData

jaspar_pwm_atsnp = list()

for (i in c(1:length(jaspar_pwm_list))){
  temp_tf = jaspar_pwm_list[[i]]
  temp_tf_name = temp_tf@name
  temp_tf_matrix = temp_tf@profileMatrix
  new_matrix = apply(temp_tf_matrix,2,function(x)(x/sum(x)))
  jaspar_pwm_atsnp[[temp_tf_name]] = t(new_matrix)
}

print("JASPAR loaded!")

output_dir = args[["output_dir"]]
#output_dir = "~/scratch/Final_pipline/Final_new_pip/results/MS/10x/"

ci_gwas_dir = args[["ci_gwas_dir"]]
#ci_gwas_dir = paste0(output_dir,"gwas_CI.txt")

ci_gwas_data = read.table(file=ci_gwas_dir, sep="\t", header=TRUE)

if(sum(c("id","chr","pos","PPA","chunk","a1","a2") %in% colnames(ci_gwas_data)) < 7){
  req_list = c("id","chr","pos","PPA","chunk","a1","a2")
  req_not_found = req_list[which(!(req_list %in% colnames(ci_gwas_data)))]
  stop(paste0("Required columns missing in CI SNPs:",req_not_found))
  
}


head(ci_gwas_data)
genome_built = args[["genome_built"]]
#genome_built = "hg19"

print("Getting CI SNPs ready for motif analysis!")

if (genome_built == "hg19"){
  library(BSgenome.Hsapiens.UCSC.hg19)
  eff_snp = ci_gwas_data
  #reg_snp_ranges = IRanges(start = eff_snp$pos, end = eff_snp$pos)
  #reg_snp_granges = GRanges(seqnames = eff_snp$chr, ranges = reg_snp_ranges)
  
  #snp_seq = as.character(getSeq(BSgenome.Hsapiens.UCSC.hg19, reg_snp_granges))
  #other_snp_seq = chartr("ATGC","TACG",snp_seq)
  
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
  #reg1_snp_data = LoadSNPData(genome.lib = "BSgenome.Hsapiens.UCSC.hg19",
  #                            snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh37",
  #                            snpids = ci_gwas_data$id,
  #                            half.window.size = 30, default.par = TRUE,
  #                            mutation = FALSE
  #)
  file.remove(snp_table_dir)
  
}else if(genome_built == "hg38"){
  library(BSgenome.Hsapiens.UCSC.hg38)
  
  eff_snp = ci_gwas_data
  #other_snp_seq = chartr("ATGC","TACG",snp_seq)
  
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
  #reg1_snp_data = LoadSNPData(genome.lib = "BSgenome.Hsapiens.UCSC.hg38",
  #                            snp.lib = "SNPlocs.Hsapiens.dbSNP144.GRCh38",
  #                            snpids = ci_gwas_data$id,
  #                            half.window.size = 30, default.par = TRUE,
  #                            mutation = FALSE
  #)
  file.remove(snp_table_dir)
}

print("Calculating Effect SNPs!")
results = ComputeMotifScore(jaspar_pwm_atsnp, reg1_snp_data, ncores = 2)


results_pval = ComputePValues(jaspar_pwm_atsnp, reg1_snp_data, results$motif.scores, ncores = 2, testing.mc=TRUE)
results_pval_val = results_pval$pval_diff
results_pval_val_cor = p.adjust(results_pval_val, method = "fdr")
significant_index = which(results_pval_val_cor < 0.1)

results_pval_sig = results_pval[significant_index,]
results_pval_val_cor = results_pval_val_cor[significant_index]

reg1_effect_snps = sort(results_pval_sig$snpid)

final_effect_snp_frame = data.frame(matrix(data = NA, 
                                           nrow = length(results_pval_sig$snpid),
                                           ncol = 9))
colnames(final_effect_snp_frame) = c("SNP","CHR","Pos","Locus","TF", "Raw_P_value", "FDR_corr_P_value",
                                     "log_like_ratio", "ppa")

for (i in c(1:length(results_pval_sig$snpid))){
  temp_id = results_pval_sig$snpid[i]
  ci_effect_index = which(ci_gwas_data$id == temp_id)
  temp_pos = ci_gwas_data$pos[ci_effect_index]
  temp_chr = ci_gwas_data$chr[ci_effect_index]
  temp_locus = ci_gwas_data$chunk[ci_effect_index]
  temp_TF = results_pval_sig$motif[i]
  temp_pval = results_pval_sig$pval_diff[i]
  temp_pval_corr = results_pval_val_cor[i]
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

#zero_index = which(final_effect_snp_frame$P_value == 0)
#final_effect_snp_frame$P_value[zero_index] = 1e-10


#final_ci_effect_dir = "~/rprojects_whole/ms_pbmc/gwas_data/bmi_chip/"
final_ci_effect_dir = paste0(output_dir, "Ci_effect_SNPs.txt")
write.table(file = final_ci_effect_dir, final_effect_snp_frame, col.names = TRUE, row.names = FALSE,
            quote = FALSE, sep = "\t")

