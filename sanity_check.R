suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(tools))

args = commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Helper function that determines what input file format the user used
# Function to read file based on extension
read_file <- function(file_path, nrows = 0) {
  ext <- file_ext(file_path)
  if (ext == "csv") {
    df = read.csv(file_path, header = TRUE,nrows=nrows)
  } else if (ext == "tsv") {
    df = read.delim(file_path, header = TRUE,nrows=nrows)
  } else {
    df = read.table(file_path, header = TRUE,nrows=nrows)
  }
  return(df)
}

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
check_path = function(path, verbose = 1){
  if (!file.exists(path)) {
    stop("Path does not exist: ", path)
  }
  else if (verbose == 1){
    message("Found ", path)
  }
}

# output_dir  
output_dir = args[["output_dir"]]
check_path(output_dir)

if (tolower(args[["finemap"]]) == "y") {
  # check if plink_binary and fgwas_src param exists
  if (nzchar(args[["plink2_dir"]])){ 
    check_path(args[["plink2_dir"]])
  } else {
    check_binary("plink2") # version checking?
  }
  if (nzchar(args[["fgwas_dir"]])){
    print(args[["fgwas_dir"]])
    check_path(args[["fgwas_dir"]])
  } else {
    check_binary("fgwas")
  }

  # sumstats
  sum_stats_dir = args[["sum_stats_dir"]]
  check_path(sum_stats_dir)
  sum_stats = read.table(sum_stats_dir, header = TRUE, nrows = 1)
  colnames(sum_stats) = tolower(colnames(sum_stats))
  req_list = c("snp","chr","pos","p")
  missing_cols = setdiff(req_list, colnames(sum_stats))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in GWAS summary statistics: ", paste(missing_cols, collapse = ", "))
  }
  # message("GWAS summary statistics columns present.")

  # leadsnps, reminder that I do not like comparing snps between datasets using RSIDs
  lead_snp_dir = args[["lead_snps_dir"]]
  check_path(lead_snp_dir)
  lead_snp = read.table(lead_snp_dir, header = TRUE, nrows = 1)
  colnames(lead_snp) = tolower(colnames(lead_snp))
  req_list = c("snp")
  missing_cols = setdiff(req_list, colnames(lead_snp))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in lead snp file: ", paste(missing_cols, collapse = ", "))
  }
  # message("Lead snp columns present.")
  
  # snp ref, reminder that the user should be able to provide their own B files which is different, so the snp_ref parameter and population parameter may not be necessary.
  # instead, they should specify the path to the B files, rather than by specifying the population
  # ex. (EAS, EUR)
  # .bed, .bim, .fam
  snp_ref_dir = args[["snp_ref_dir"]]
  prefix_name = args[["population"]]
  extension_list = c(".bed",".bim",".fam") # plink2 also have an alternative extension name for binary files, see pgen
  for (ext in extension_list){
    filename = paste0(prefix_name, ext)
    check_path(paste0(output_dir,filename))
  }
} else {
  # ci_gwas_dir
  ci_gwas_dir = args[["ci_gwas_dir"]]
  ci_gwas = read.table(ci_gwas_dir, header = TRUE, nrows = 0)
  colnames(ci_gwas) = tolower(colnames(ci_gwas))
  req_list = c("id","chr","pos","ppa","chunk")
  missing_cols = setdiff(req_list, colnames(ci_gwas))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in credible interval SNPs: ", paste(missing_cols, collapse = ", "))
  }
 # message("Credible interval gwas columns present.")
}

# genome_build
genome_build = args[["genome_build"]]
genome_build = tolower(genome_build)
req_builds = c("hg19","hg38")
if (!(genome_build %in% req_builds)) {
  stop("Invalid genome build: '", genome_build, 
       "'. Allowed values are: ", paste(req_builds, collapse = ", "))
}

# gencode_dir
gene_annot_dir = args[["gencode_dir"]]
check_path(gene_annot_dir)

if (tolower(args[["histone_mark_analysis"]]) == "y") {
  hist_mark_instruct_dir = args[["hist_mark_instruct_dir"]]
  check_path(hist_mark_instruct_dir)
  hist_mark_instruct = read.table(hist_mark_instruct_dir, nrows = -1)
  colnames(hist_mark_instruct) = tolower(colnames(hist_mark_instruct))
  req_list = c("chromhmm_dir","atac_cell_types","chromhmm_cell_types")
  missing_cols = setdiff(req_list, colnames(hist_mark_instruct))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in histone mark instructions file: ", paste(missing_cols, collapse = ", "))
  }
  # check if contents of histone marker instruct directories exist
  for (i in hist_mark_instruct$chromhmm_dir){
    check_path(i, verbose = 0)
  }
  
  # message("Histone mark instructions columns present.")
}

if (tolower(args[["hic_analysis"]]) == "y") {
  hic_instruct_dir = args[["hic_instruct_dir"]]
  check_path(hic_instruct_dir)
  hic_instruct = read_file(hic_instruct_dir, nrows = -1)
  colnames(hic_instruct) = tolower(colnames(hic_instruct))
  req_list = c("hic_dir","genes_present","bulk","atac_cell_types","hic_cell_types")
  missing_cols = setdiff(req_list, colnames(hic_instruct))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in HI-C instructions file: ", paste(missing_cols, collapse = ", "))
  }
  # check if contents of hic instruct directories exist
  for (i in hic_instruct$hic_dir){
    check_path(i, verbose = 0)
  }
  # message("HI-C instructions columns present.")
}

# need to check if cell type columns matches
if (tolower(args[["eqtl_analysis"]]) == "y") {
  eqtl_instruct_dir = args[["eqtl_instruct_dir"]]
  eqtl_instruct_dir = "/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/eqtl_instructions.tsv"
  check_path(eqtl_instruct_dir)
  eqtl_instruct = read_file(eqtl_instruct_dir, nrows= -1)
  colnames(eqtl_instruct) = tolower(colnames(eqtl_instruct))
  req_list = c("eqtl_dir","atac_cell_types")
  missing_cols = setdiff(req_list, colnames(eqtl_instruct))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in eQTL instructions file: ", paste(missing_cols, collapse = ", "))
  }

  # check if contents of eqtl instruct directories exist
  for (i in eqtl_instruct$eqtl_dir){
    check_path(i, verbose = 0)
  }
  # message("eQTL instructions columns present.")
}

tf_setting = args[["tf_expr_analysis"]]
tf_setting = tolower(tf_setting)
tf_list = c("atac","rna","both","none")
if (!(tf_setting %in% tf_list)) {
  stop("Invalid TF option: '", tf_setting, 
       "'. Allowed values are: ", paste(tf_list, collapse = ", "))
}

if (tolower(args[["tf_expr_analysis"]]) %in% c("rna", "both")) {
  scrna_instruct_dir = args[["scrna_instruct_dir"]]
  check_path(scrna_instruct_dir)
  scrna_instruct = read_file(scrna_instruct_dir, nrows=-1)
  colnames(scrna_instruct) = tolower(colnames(scrna_instruct))
  req_list = c("scrna_dir","atac_cell_types","rna_cell_types","one_cell_type_seurat","matrix_loc")
  missing_cols = setdiff(req_list, colnames(scrna_instruct))
  if (length(missing_cols) > 0) {
    stop("Required columns missing in scrna instructions file: ", paste(missing_cols, collapse = ", "))
  }
  
  # check if contents of eqtl instruct directories exist
  for (i in scrna_instruct$scrna_dir){
    check_path(i, verbose = 0)
  }
  # message("scrna instructions columns present.")
}

# Settings output
settings_msg <- paste(
  "------------------------------------\n",
  "| RegSCOUT Settings                |\n",
  "|                                  |\n",
  "| Genome Build            : %-6s |\n",
  "|                                  |\n",
  "| Finemapping             : %-6s |\n",
  "| TF Expression Analysis  : %-6s |\n",
  "| Histone Marker Analysis : %-6s |\n",
  "| Hi-C Analysis           : %-6s |\n",
  "| eQTL Analysis           : %-6s |\n",
  "------------------------------------\n",
  sep = ""
)

cat(sprintf(settings_msg,
  tolower(args[["genome_build"]]),
  if (tolower(args[["finemap"]]) == "y") toupper(args[["finemap"]]) else "N",
  tolower(args[["tf_expr_analysis"]]),
  toupper(args[["histone_mark_analysis"]]),
  toupper(args[["hic_analysis"]]),
  toupper(args[["eqtl_analysis"]])
))