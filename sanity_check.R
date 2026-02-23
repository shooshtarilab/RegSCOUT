suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(xlsx))

# check for libraries

args = commandArgs(trailingOnly = TRUE, asValues = TRUE)
# write.xlsx(test, "/workspaces/RegSCOUT/RA/inputs/instruction_files/hic_instructions.xlsx", row.names = F)
# Helper function that determines what input file format the user used
# Function to read file based on extension

read_file <- function(file_path, nrows = 0) { # can process xlsx, csv, tsv
  ext <- file_ext(file_path)
  if (ext == "xlsx") {
    if (nrows == Inf) nrows = NULL
    df = read.xlsx(file_path, rowIndex = nrows, sheetIndex = 1)
  } else {
    df = fread(file_path, nrows=nrows) # should handle any kind of separator
  }
  # xlsx
  return(df)
}

check_binary = function(binary) {
  path = Sys.which(binary)
  if (nchar(path) == 0) {
    stop(
      sprintf(
        "Binary '%s' not found in PATH. Please add it to PATH or specify path manually using --plink2_bin or --fgwas_src.",
        binary
      ),
      call. = F
    )
  } else {
   #  message(sprintf("Found %s at '%s'", binary, path))
  }
}

# also add check if directory is read/writable?
check_path = function(path, name, dir=F, verbose = 0){
  if (dir){
    if (!dir.exists(path)) {
      stop("Directory does not exist '", path,"'", call. = FALSE)
    }   
  } else {
    if (!file.exists(path)) {
      stop("File does not exist '", path,"'", call. = FALSE)
    }
  }
  if (verbose == 1){
    message(sprintf("%s found at'%s'", name, path))
  }
}

check_col_diff = function(file_path, required_cols, name){
  header = read_file(file_path, nrows= 0)
  missing_cols = setdiff(tolower(required_cols), tolower(colnames(header)))
  is_missing = !(tolower(required_cols) %in% tolower(colnames(header)))
  if (any(is_missing)) {
      stop("Error: Missing columns: ", paste(required_cols[is_missing], collapse = ", "))
  }
}

check_instruction_contents = function(instruct_f, dirname){
  instruction_file = read_file(instruct_f, Inf)
  dir_list = instruction_file[[dirname]]
  lapply(dir_list, check_path, verbose = 0)

  if (dirname  == "eqtl_dir"){
    dir_list = instruction_file[instruction_file[["tabix"]] == TRUE, ][["eqtl_dir"]]
    dir_list = lapply(dir_list, paste0, ".tbi")
    lapply(dir_list, check_path, verbose = 0)
  }
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
  "| Histone Mark Analysis   : %-6s |\n",
  "| Hi-C Analysis           : %-6s |\n",
  "| eQTL Analysis           : %-6s |\n",
  "------------------------------------\n",
  sep = ""
)
cat(do.call(sprintf, c(
  list(settings_msg, tolower(args[["genome_build"]])),
  lapply(args[c("finemap", "tf_expr_analysis", "histone_mark_analysis", "hic_analysis", "eqtl_analysis")], 
         \(x) if(identical(tolower(x), "y")) "Y" else "N")
)))


message("\n --- Checking for output directory ---")

output_dir = args[["output_dir"]]
check_path(output_dir, "Output directory", T, verbose=1)

# create folder structure
figure_dir = paste0(output_dir,"/figures")
if (dir.exists(figure_dir)) {
  unlink(figure_dir, recursive = TRUE, force = TRUE)
}
dir.create(figure_dir)
message("Pass")

# Check finemapping requirements
if (tolower(args[["finemap"]]) == "y") {
  # Check binaries
  message("\n--- Checking finemapping requirements ----")
  # check if plink_binary and fgwas_src param exists
  if (nzchar(args[["plink2_dir"]])){ 
    check_path(args[["plink2_dir"]])
  } else {
    check_binary("plink2")
  }
  if (nzchar(args[["fgwas_dir"]])){
    check_path(args[["fgwas_dir"]])
  } else {
    check_binary("fgwas")
  }

  # Check finemapping files
  check_path(args[["sum_stats_file"]], name = "GWAS summary statistics file")
  check_col_diff(args[["sum_stats_file"]], c("chr","p","pos","snp"), "GWAS summary statistics file")
  check_path(args[["lead_snps_file"]], name = "Lead SNPs file") 
  check_col_diff(args[["lead_snps_file"]], c("snp"), "Lead SNPs file")

  check_path(args[["snp_ref_dir"]], name="SNP reference file directory", dir=T)
  snp_ref_file = list.files(path= args[["snp_ref_dir"]], pattern = "\\.(bed|fam|bim)", ignore.case = T)
  required_files = paste0(args[["population"]], c(".bim",".fam",".bed"))
  missing <- setdiff(required_files, snp_ref_file)

  if (length(missing) > 0) {
    stop(paste("Error: Missing required PLINK files:", paste(missing, collapse = ", ")), call. = FALSE)
  } 
  message("Pass")
} else {
    message("\n--- Checking finemapping data ----")
    check_path(args[["ci_gwas_file"]])
    check_col_diff(args[["ci_gwas_file"]],c("id","pos","chr","a1","a2","chunk","ppa"), name= "Finemapped SNP file")
    
    if (nzchar(args[["loci_info_file"]])){
      check_path(args[["loci_info_file"]])
      check_col_diff(args[["loci_info_file"]],c("chunk"),name= "Loci information file")

      # Warning for "locus_chr", "locus_start", "locus_end"
      header = read_file(args[["loci_info_file"]], nrows= 0)
      loci_cols = c("locus_chr", "locus_start", "locus_end")
      missing_cols = setdiff(tolower(loci_cols), tolower(colnames(header)))
      is_missing = !(tolower(loci_cols) %in% tolower(colnames(header)))
      if (any(is_missing)) {
        warning("Loci info file is missing: ", paste(loci_cols[is_missing], collapse = ", "), ". This information will not be included in final tables.")
      }
    }
  message("Pass")
}

# Check finemapping date
# Check mode peak_table
if (tolower(args[["mode"]] == "peak_table")){
  message("\n--- Checking mode setting peak_table ---")
  check_path(args[["peak_cell_file"]])
  check_col_diff(args[["peak_cell_file"]], c("chr","start", "end", "cell_types"), name = "Peak-Cell file") 

  check_path(args[["cicero_dir"]], dir=T)
  cicero_files = list.files(path= args[["cicero_dir"]], pattern = "\\.filtered_coaccessible_sites4s.txt", ignore.case = T, full.names = T)
  cell_name = basename(cicero_files) |> sub("\\..*","", x = _)
  message(sprintf("Cicero cell types found: %s", paste(cell_name, collapse = ", ")))
  for (file in cicero_files){
    cell_name = basename(file) |> sub("\\..*","", x = _)
    check_col_diff(file, c("Peak1","Peak2", "coaccess"), name = paste0("Cicero file ", cell_name)) 
  }

  genome_build = tolower(args[["genome_build"]])
  req_builds = c("hg19","hg38")
  if (!(genome_build %in% req_builds)) {
    stop("Invalid genome build: '", genome_build, 
         "'. Allowed values are: ", paste(req_builds, collapse = ", "))
  }
  check_path(args[["gencode_file"]])
  gencode_file = tolower(basename(args[["gencode_file"]]))

  if (endsWith(gencode_file, ".gff3.gz") || endsWith(gencode_file, ".gff3" || endsWith(gencode_file, ".gtf.gz") || endsWith(gencode_file, ".gtf"))){
  } else {
    stop(paste("Error: Expected GENCODE file to end with gff3, gff3.gz, gtf or gtf.gz"), call. = FALSE)
  }

  if (!nzchar(args[["jaspar_mtx_file"]])){
    warning("--jaspar_mtx_file is unset. Using pipeline default file.")
  }
  message("Pass")
}

# Check mode atac_obj
if (tolower(args[["mode"]] == "atac_obj")){
  message("\n--- Checking mode setting atac_obj---")
  check_path(args[["seurat_obj_file"]])
  if (tolower(file_ext(args[["seurat_obj_file"]])) == "rds"){
    seurat_obj <- readRDS(args[["seurat_obj_file"]])
  } else if (tolower(file_ext(args[["seurat_obj_file"]])) == "rdata"){
    seurat_obj = load(args[["seurat_obj_file"]])
    seurat_obj= get(seurat_obj)
  }
  message("Pass")
}

  # Check for assay named peaks
  # from counts matrix
  # cell type labels, acessible using Idents

# Checking instruction file directory
if (tolower(args[["histone_mark_analysis"]]) == "y" || tolower(args[["hic_analysis"]]) == "y" || tolower(args[["eqtl_analysis"]]) == "y"|| tolower(args[["tf_expr_analysis"]]) %in% c("atac", "both", "rna")){
  histone_f= hic_f= eqtl_f= scrna_f = F
  message("\n--- Checking instruction file directory ----") 
  instruct_dir = args[["instruction_file_dir"]]
  if (!nzchar(instruct_dir)){
    stop("--instruction_file_dir is unset.", call. = FALSE)
  }
  check_path(instruct_dir, dir=T)
  instruct_files = list.files(path = instruct_dir, pattern = "\\.(csv|tsv|txt)", all.files= F, full.names = T, ignore.case = T)
  file_header_lookup = c(
    "atac_cell_types,cell_sorted,genes_present,hic_cell_types,hic_dir,signif_threshold" = "hic_instruct",
    "atac_cell_types,eqtl_dir,tabix" = "eqtl_instruct",
    "atac_cell_types,chromhmm_cell_types,chromhmm_dir" = "histone_mark_instruct",
    "atac_cell_types,matrix_loc,rna_cell_types,scrna_dir,one_cell_type_seurat" = "scrna_instruct"
  )

  bash_vars <- c()

  for (f in instruct_files){
    f_basename = basename(f)
    header = try(fread(f, nrows = 0), silent = T)
    if(inherits(header, "try-error")){
      warning(sprintf("%s cannot be read.",f))
      next
    }

    file_cols_vec <- colnames(header) |> trimws() |> tolower()
    match_idx <- sapply(names(file_header_lookup), function(lookup_str) {
      lookup_cols_vec <- strsplit(lookup_str, ",")[[1]]
      all(lookup_cols_vec %in% file_cols_vec) # return T if every required col from the lookup entry is present in file.
    })


    if (any(match_idx)) {
      var_name <- file_header_lookup[[which(match_idx)[1]]] 
      switch(var_name,
        "histone_mark_instruct" = {
          if (!nzchar(args[["histone_mark_analysis"]]) ||  tolower(args[["histone_mark_analysis"]]) != "y") {
            #warning("Histone marker instruction file detected but analysis is not enabled.")
            next
          } else {
            message(sprintf("Found histone marker instruction file '%s'", f_basename))  
            check_instruction_contents(f, "chromHMM_dir")
            histone_f = T
          }
        },
        "scrna_instruct" = {
          if (!nzchar(args[["tf_expr_analysis"]]) || ! tolower(args[["tf_expr_analysis"]]) %in% c("both", "rna")) {
            next
          } else {
            message(sprintf("Found scRNA instruction file '%s'", f_basename))  
            check_instruction_contents(f, "scrna_dir")
            scrna_f = T
          }
        },
        "hic_instruct" = {
          if (!nzchar(args[["hic_analysis"]]) ||  tolower(args[["hic_analysis"]]) != "y") {
            #warning("Hi-C instruction file detected but analysis is not enabled.")
            next
          } else {
            message(sprintf("Found Hi-C instruction file '%s'", f_basename))  
            check_instruction_contents(f, "hic_dir")
            hic_f = T
          }
        },
        "eqtl_instruct" = {
          if (!nzchar(args[["eqtl_analysis"]]) || tolower(args[["eqtl_analysis"]]) != "y") {
            # warning("eQTL instruction file detected but analysis is not enabled.")
            next
          } else {
            message(sprintf("Found eQTL instruction file '%s'",f_basename)) 
            check_instruction_contents(f, "eqtl_dir")
            eqtl_f = T
          }
        }
      )
      bash_vars <- c(bash_vars, paste0("export ", var_name, "='", f, "'"))
    }
  }

  if (args[["hic_analysis"]] == "y" && !hic_f) {
    stop("Hi-C analysis enabled but instruction file not found. Check column headers.", call. = FALSE)
  }
    if (args[["eqtl_analysis"]] == "y" && !eqtl_f) {
    stop("eQTL analysis enabled but instruction file not found. Check column headers.", call. = FALSE)
  }
    if (args[["histone_mark_analysis"]] == "y" && !histone_f) {
    stop("Hi-C analysis enabled but instruction file not found. Check column headers.", call. = FALSE)
  }
    if (args[["tf_expr_analysis"]] %in% c("both","rna") && !scrna_f) {
    stop("Hi-C analysis enabled but instruction file not found. Check column headers.", call. = FALSE)
  }

  message("Pass")
  writeLines(bash_vars, "temp_instruct_vars.sh")
}