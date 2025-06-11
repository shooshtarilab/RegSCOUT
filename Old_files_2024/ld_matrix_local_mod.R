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
