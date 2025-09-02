suppressPackageStartupMessages(library(R.utils))

args = commandArgs(trailingOnly = TRUE, asValues = TRUE)

check_binary = function(binary) {
  path = Sys.which(binary)
  
  if (nchar(path) == 0) {
    stop(
      sprintf(
        "Required binary '%s' not found in PATH. Please install it or add it to PATH.",
        binary
      ),
      call. = FALSE
    )
  } else {
    message(sprintf("Found '%s' at: %s", binary, path))
    return(invisible(TRUE))
  }
}

check_path = function(path){
  if (!file.exists(path)) {
    stop("Directory does not exist: ", path)
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

  

}




# check if directory is read/writable?

# output_dir
output_dir = args[["output_dir"]]
check_path(output_dir)


# check ci_gwas_dir (finemapped snps)
# A tab-separated file which contains the results of fine-mapping and filtering for CI_thr.
# SNPs on sex chromosomes and with PPA values <= 0.01 have also been filtered out. why?
# These are the credible interval SNPs.

# ci_gwas_dir
ci_gwas_dir = args[["ci_gwas_dir"]]
ci_gwas = read.table(ci_gwas_dir, header = TRUE, nrows = 0)
colnames(ci_gwas) = tolower(colnames(ci_gwas))
req_list = c("id","chr","pos","ppa","chunk")
missing_cols = setdiff(req_list, colnames(ci_gwas))
if (length(missing_cols) > 0) {
  stop("Required columns missing in CI SNPs: ", paste(missing_cols, collapse = ", "))
}
message("CI_gwas columns available.")


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

