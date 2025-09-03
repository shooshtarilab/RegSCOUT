suppressPackageStartupMessages(library(R.utils))

args = commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Helper function that determines what input file format the user used
# Supports csv, tsv, txt
library(tools)

# Function to read file based on extension
read_file <- function(file_path) {
  ext <- file_ext(file_path)
  
  # Decide based on extension
  if (ext == "csv") {
    df = read.csv(file_path, header = TRUE)
  } else if (ext == "tsv") {
    df = read.delim(file_path, header = TRUE)
  } else if (ext == "txt") {
    df = read.table(file_path, sep = " ", header = TRUE)
  } else {
    stop("Unsupported file format: '.", ext, "' in ", file_path)
  }
  return(df)
}

# Example usage
file_path <- "/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hic_instructions.tsdv"   # user specifies the file path
df <- read_file(file_path)
print(df)



check_binary = function(binary) {
  path = Sys.which(binary)
  if (nchar(path) == 0) {
    stop(
      sprintf(
        "Required binary '%s' not found in PATH. Please add add it to PATH or specify path manually using --plink2_bin or --fgwas_src.",
        binary
      ),
      call. = FALSE
    )
  } else {
    message(sprintf("Found '%s' at: %s", binary, path))
  }
}

# also add check if directory is read/writable?
check_path = function(path){
  if (!file.exists(path)) {
    stop("Path does not exist: ", path)
  }
  else{
    message("Found ", path)
  }
}

if (tolower(args[["finemap"]]) == "y") {
  # check if plink_binary and fgwas_src param exists
  if (!is.null(args[["plink2_bin"]])){
    check_path(args[["plink2_bin"]])
  } else {
    check_binary("plink2") # version checking?
  }
  if (!is.null(args[["fgwas_src"]])){
    check_path(args[["fgwas_src"]])
  } else {
    check_binary("fgwas")
  }

  # sumstats
  sum_stats_dir = args[["sum_stats"]]
  check_path(sum_stats_dir)
  sum_stats = read.table(sum_stats_dir, header = TRUE, nrows = 1)
  colnames(sum_stats) = tolower(colnames(sum_stats))
  req_list = c("snp","chr","pos","p")
  missing_cols = setdiff(req_list, colnames(sum_stats))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in GWAS summary statistics: ", paste(missing_cols, collapse = ", "))
  }
  message("GWAS summary statistics columns present.")

  # leadsnps, reminder that I do not like comparing snps between datasets using RSIDs
  lead_snp_dir = args[["lead_snp"]]
  check_path(lead_snp_dir)
  lead_snp = read.table(lead_snp_dir, header = TRUE, nrows = 1)
  colnames(lead_snp) = tolower(colnames(lead_snp))
  req_list = c("snp")
  missing_cols = setdiff(req_list, colnames(lead_snp))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in lead snp file: ", paste(missing_cols, collapse = ", "))
  }
  message("Lead snp columns present.")

  # snp ref, reminder that the user should be able to provide their own B files which is different, so the snp_ref parameter and population parameter may not be necessary.
  # instead, they should specify the path to the B files, rather than by specifying the population
  # ex. (EAS, EUR)
  # .bed, .bim, .fam
  # snp_ref_dir = args[["SNP_ref"]]
  # check_path(snp_ref_dir)
  # snp_ref = read.table(snp_ref_dir, header = TRUE, nrows = 1)
  # colnames(snp_ref) = tolower(colnames(snp_ref))
  # req_list = c("EUR, EAS, AFR") # what else?
  #  missing_cols = setdiff(req_list, colnames(lead_snp))
  # if (length(missing_cols) > 0) {
  #   stop("Required columns missing in lead snp file: ", paste(missing_cols, collapse = ", "))
  # }
  # message("Lead snp columns present.")
}

# output_dir
output_dir = args[["output_dir"]]
check_path(output_dir)

# ci_gwas_dir
ci_gwas_dir = args[["ci_gwas_dir"]]
ci_gwas = read.table(ci_gwas_dir, header = TRUE, nrows = 0)
colnames(ci_gwas) = tolower(colnames(ci_gwas))
req_list = c("id","chr","pos","ppa","chunk")
missing_cols = setdiff(req_list, colnames(ci_gwas))
if (length(missing_cols) > 0) {
  stop("Required columns missing in credible interval SNPs: ", paste(missing_cols, collapse = ", "))
}
message("Credible interval gwas columns present.")


# genome_built 
genome_build = args[["genome_built"]]
genome_build = tolower(genome_build)
req_builds = c("hg19","hg38")
if (!(genome_build %in% req_builds)) {
  stop("Invalid genome_build: '", genome_build, 
       "'. Allowed values are: ", paste(req_builds, collapse = ", "))
}
message("Using genome_build: ", genome_build)

# genecode_dir
gene_annot_dir = args[["gencode_dir"]]
check_path(gene_annot_dir)

if (tolower(args[["histone_mark_analysis"]]) == "y") {
  hist_mark_instruct_dir = args[["hist_mark_instruct_dir"]]
  check_path(hist_mark_instruct_dir)
  hist_mark_instruct = read.table(hist_mark_instruct_dir, header = TRUE, nrows = 0)
  colnames(hist_mark_instruct) = tolower(colnames(hist_mark_instruct))
  print(head(hist_mark_instruct))
  req_list = c("chromhmm_dir","atac_cell_types","chromhmm_cell_types")
  missing_cols = setdiff(req_list, colnames(hist_mark_instruct))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in credible interval SNPs: ", paste(missing_cols, collapse = ", "))
  }
  message("Histone mark instructions columns present.")
}

if (tolower(args[["histone_mark_analysis"]]) == "y") {
  hist_mark_instruct_dir = args[["hist_mark_instruct_dir"]]
  check_path(hist_mark_instruct_dir)
  hist_mark_instruct = read.table(hist_mark_instruct_dir, header = TRUE, nrows = 0)
  print(hist_mark_instruct)
  colnames(hist_mark_instruct) = tolower(colnames(hist_mark_instruct))
  req_list = c("chromHMM_dir","atac_cell_types","chromHMM_cell_types")
  missing_cols = setdiff(req_list, colnames(hist_mark_instruct))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in credible interval SNPs: ", paste(missing_cols, collapse = ", "))
  }
  message("Histone mark instructions columns present.")
}

output_dir = "/home/ubunkun/Lab/RA_project/RegSCOUT/EUR"
