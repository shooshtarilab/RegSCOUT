library(R.utils)

args <- commandArgs(trailingOnly = TRUE, asValues = TRUE)

output_dir = args[["output_dir"]]
#output_dir = "/home/nader95/scratch/Final_pipline/Final_new_pip/results/MS/10x/"
fgwas_src = args[["fgwas_src"]]
#fgwas_src = "/home/nader95/scratch/Final_pipline/Final_new_pip/code/fgwas-0.3.6/src/fgwas"
fgwas_file = paste0(output_dir,"final_gwas_data.txt")
ci_suff = "CI"
ci_files = paste0(output_dir,ci_suff)
print(fgwas_file)
shell <- ifelse(Sys.info()['sysname'] == "Windows", "cmd", "sh")
fun <- paste0(
    shQuote(fgwas_src, type=shell),
    " -i ", shQuote(fgwas_file, type=shell),
    " -o ", shQuote(ci_files, type=shell),
    " -fine ",
    " -print "
)
system(fun)

ci_file_bfs = paste0(ci_files,".bfs.gz")

fun <- paste0(
    shQuote("gzip", type=shell),
    " -d ",
    shQuote(ci_file_bfs, type=shell)
)

system(fun)

ci_files_bfs = paste0(ci_files,".bfs")
ci_data = read.table(ci_files_bfs, sep = " ", header = TRUE)
fgwas_data = read.table(fgwas_file, sep = "\t", header = TRUE)

ci_data[["a1"]] = fgwas_data$A1
ci_data[["a2"]] = fgwas_data$A2

print("Extracting CI SNPs")
ci_th = as.numeric(args[["CI_thr"]])
#ci_th = 0.95
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

chr_list = paste0("chr", c(1:22))
ci_gwas_data = ci_gwas_data[ci_gwas_data$chr %in% chr_list,]

ci_gwas_data = ci_gwas_data[ci_gwas_data$PPA > 0.01,]

ci_dir = paste0(output_dir,"gwas_CI.txt")
write.table(ci_gwas_data, file = ci_dir, col.names = TRUE, sep="\t",
            row.names = FALSE, quote = FALSE)
