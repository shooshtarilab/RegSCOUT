suppressPackageStartupMessages(library(R.utils))
suppressPackageStartupMessages(library(tools))
suppressPackageStartupMessages(library(data.table))
# check for libraries

args = commandArgs(trailingOnly = TRUE, asValues = TRUE)

# Helper function that determines what input file format the user used
# Function to read file based on extension
read_file <- function(file_path, nrows = 1) {
  ext <- file_ext(file_path)
  if (ext == "csv") {
    df = read.csv(file_path, nrows=nrows)
  } else if (ext == "tsv") {
    df = fread(file_path, nrows=nrows)
  } else {
    df = fread(file_path, nrows=nrows)
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
}

# Check finemapping data
if (tolower(args[["finemap"]]) == "n" || !nzchar(args[["finemap"]])) {
  message("\n--- Checking finemapping data ----")
  check_path(args[["ci_gwas_file"]])
  check_col_diff(args[["ci_gwas_file"]],c("id","pos","chr","a1","a2","chunk","ppa"), name= "Finemapped SNP file")
  message("Pass")
}

# Check mode peak_table
if (tolower(args[["mode"]] == "peak_table")){
  message("\n--- Checking mode setting peak_table---")
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

    check_path(args[["gencode_file"]])
    gencode_file = tolower(basename(args[["gencode_file"]]))

    if (endsWith(gencode_file, ".gff3.gz") || endsWith(gencode_file, ".gff3" || endsWith(gencode_file, ".gtf.gz") || endsWith(gencode_file, ".gtf"))){
    } else {
      stop(paste("Error: Expected GENCODE file to end with gff3, gff3.gz, gtf or gtf.gz"), call. = FALSE)
    }

    if (!nzchar(args[["jaspar_mtx_file"]])){
      warning("--jaspar_mtx_file is unset. Using pipeline default file.")
    }
    
  }
  message("Pass")
}

# Check mode atac_obj
if (tolower(args[["mode"]] == "atac_obj")){
  message("\n--- Checking mode setting atac_obj---")
  if (tolower(file_ext(args[["seurat_obj_file"]])) == "rds"){
    seurat_obj <- readRDS(args[["seurat_obj_file"]])
  } else if (tolower(file_ext(args[["seurat_obj_file"]])) == "rdata"){
    seurat_obj = load(args[["seurat_obj_file"]])
    seurat_obj= get(seurat_obj)
  }
}

  # Check for assay named peaks
  # from counts matrix
  # cell type labels, acessible using Idents

  genome_build = tolower(args[["genome_build"]])
  req_builds = c("hg19","hg38")
  if (!(genome_build %in% req_builds)) {
    stop("Invalid genome build: '", genome_build, 
         "'. Allowed values are: ", paste(req_builds, collapse = ", "))
  
  check_path(args[["gencode_file"]])
  gencode_file = tolower(basename(args[["gencode_file"]]))

  if (endsWith(gencode_file, ".gff3.gz") || endsWith(gencode_file, ".gff3" || endsWith(gencode_file, ".gtf.gz") || endsWith(gencode_file, ".gtf"))){
  } else {
     stop("Error: Expected GENCODE file to end with gff3, gff3.gz, gtf or gtf.gz", call. = FALSE)
  }

  if (!nzchar(args[["jaspar_mtx_file"]])){
    warning("--jaspar_mtx_file is unset. Using pipeline default file.")
  }
  message("Pass")
}

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
            histone_f = T
          }
        },
        "scrna_instruct" = {
          if (!nzchar(args[["tf_expr_analysis"]]) || ! tolower(args[["tf_expr_analysis"]]) %in% c("both", "rna")) {
            next
          } else {
            message(sprintf("Found scRNA instruction file '%s'", f_basename))  
            scrna_f = T
          }
        },
        "hic_instruct" = {
          if (!nzchar(args[["hic_analysis"]]) ||  tolower(args[["hic_analysis"]]) != "y") {
            #warning("Hi-C instruction file detected but analysis is not enabled.")
            next
          } else {
            message(sprintf("Found Hi-C instruction file '%s'", f_basename))  
            hic_f = T
          }
        },
        "eqtl_instruct" = {
          if (!nzchar(args[["eqtl_analysis"]]) || tolower(args[["eqtl_analysis"]]) != "y") {
            # warning("eQTL instruction file detected but analysis is not enabled.")
            next
          } else {
            message(sprintf("Found eQTL instruction file '%s'",f_basename)) 

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
}
message("Pass")
writeLines(bash_vars, "temp_instruct_vars.sh")
  
#   # check if contents of eqtl instruct directories exist
#   for (i in 1:nrow(eqtl_instruct)){
#     check_path(eqtl_instruct$eqtl_dir[i], verbose = 0)
#     # check if index file exists if tabix option is T
#     if ("tabix" %in% colnames(eqtl_instruct) && as.logical(toupper(eqtl_instruct$tabix[i]))) {
#       check_path(paste0(eqtl_instruct$eqtl_dir[i],".tbi"), verbose = 0)
#     } else { # for non tabix file, ensure colnames exist
#       # user must either provide "snp" OR "chr+pos" OR both
#       eqtl_data = read_file(eqtl_instruct$eqtl_dir[i])
#       colnames(eqtl_data) = tolower(colnames(eqtl_data))
#       req_snp <- "snp"
#       req_chrpos <- c("chr", "pos")
#       has_snp <- req_snp %in% colnames(eqtl_data)
#       has_chrpos <- all(req_chrpos %in% colnames(eqtl_data))
#       if (!(has_snp || has_chrpos)) {
#         stop('Must have either "chr" & "pos" OR "snp" columns (or both) in: ', basename(eqtl_instruct$eqtl_dir[i]))
#       } else if (has_snp && !has_chrpos) { # warn user that they only have a snp column
#         warning('Only "snp" column found (no "chr" & "pos") in: ', basename(eqtl_instruct$eqtl_dir[i]))
#       }
#     }
#   }
#   # message("eQTL instructions columns present.")
# }

# # lead snps file only needs snp??? + whatever columns

# # finemap file: id, pos, chr, chunk, ppa, a1, a2, optional: lead_snp, locus_chr, locus_start, locus_end

# # must have one reference and alternative allele, can include it in a warning

# # cicero files + cell_peak only for peak-table mode 'filtered_coaccessible_sites4s.txt'
# # cicero files are output by a rscript
# # mode atac_obj
# # loci info?




# # -------------------
# head(gwas_sum_stats)
#   # sumstats
#   sum_stats_dir = args[["sum_stats_dir"]]
#   check_path(sum_stats_dir)
#   sum_stats = read.table(sum_stats_dir, header = TRUE, nrows = 1)
#   colnames(sum_stats) = tolower(colnames(sum_stats))
#   req_list = c("snp","chr","pos","p")
#   missing_cols = setdiff(req_list, colnames(sum_stats))
#   if (length(missing_cols) > 0) {
#     stop("Required columns missing in GWAS summary statistics: ", paste(missing_cols, collapse = ", "))
#   }
#   # message("GWAS summary statistics columns present.")

#   # leadsnps, reminder that I do not like comparing snps between datasets using RSIDs
#   lead_snp_dir = args[["lead_snps_dir"]]
#   check_path(lead_snp_dir)
#   lead_snp = read.table(lead_snp_dir, header = TRUE, nrows = 1)
#   colnames(lead_snp) = tolower(colnames(lead_snp))
#   req_list = c("snp")
#   missing_cols = setdiff(req_list, colnames(lead_snp))
#   if (length(missing_cols) > 0) {
#     stop("Required columns missing in lead snp file: ", paste(missing_cols, collapse = ", "))
#   }
#   # message("Lead snp columns present.")
  
#   # snp ref, reminder that the user should be able to provide their own B files which is different, so the snp_ref parameter and population parameter may not be necessary.
#   # instead, they should specify the path to the B files, rather than by specifying the population
#   # ex. (EAS, EUR)
#   # .bed, .bim, .fam
#   snp_ref_dir = args[["snp_ref_dir"]]
#   prefix_name = args[["population"]]
#   extension_list = c(".bed",".bim",".fam") # plink2 also have an alternative extension name for binary files, see pgen
#   for (ext in extension_list){
#     filename = paste0(prefix_name, ext)
#     check_path(paste0(snp_ref_dir,filename))
#   }
# } else {
#   # ci_gwas_dir
#   ci_gwas_dir = args[["ci_gwas_dir"]]
#   ci_gwas = read.table(ci_gwas_dir, header = TRUE, nrows = 0)
#   colnames(ci_gwas) = tolower(colnames(ci_gwas))
#   req_list = c("id","chr","pos","ppa","chunk")
#   missing_cols = setdiff(req_list, colnames(ci_gwas))
#   if (length(missing_cols) > 0) {
#     stop("Required columns missing in credible interval SNPs: ", paste(missing_cols, collapse = ", "))
#   }
#   if (nrow(ci_gwas) == 0){
#     stop("Provided finemapped data (ci_gwas_dir) have no results.")
#   }
#  # message("Credible interval gwas columns present.")
# }

# if (nzchar(args[["mode"]]) && tolower(args[["mode"]]) %in% c("atac_obj", "peak_table")) {
#   # genome_build
#   genome_build = args[["genome_build"]]
#   genome_build = tolower(genome_build)
#   req_builds = c("hg19","hg38")
#   if (!(genome_build %in% req_builds)) {
#     stop("Invalid genome build: '", genome_build, 
#          "'. Allowed values are: ", paste(req_builds, collapse = ", "))
#   }
  
#   # gencode_dir
#   gene_annot_dir = args[["gencode_dir"]]
#   check_path(gene_annot_dir)
  
#   if (tolower(args[["histone_mark_analysis"]]) == "y") {
#     hist_mark_instruct_dir = args[["hist_mark_instruct_dir"]]
#     check_path(hist_mark_instruct_dir)
#     hist_mark_instruct = read.table(hist_mark_instruct_dir, header = T, nrows = -1)
#     colnames(hist_mark_instruct) = tolower(colnames(hist_mark_instruct))
#     req_list = c("chromhmm_dir","atac_cell_types","chromhmm_cell_types")
#     missing_cols = setdiff(req_list, colnames(hist_mark_instruct))
#     if (length(missing_cols) > 0) {
#       stop("Required columns missing in histone mark instructions file: ", paste(missing_cols, collapse = ", "))
#     }
#     if (nrow(hist_mark_instruct)==0){
#       stop("Histone mark instruction file is empty.")
#     }
#     # check if contents of histone marker instruct directories exist
#     for (i in hist_mark_instruct$chromhmm_dir){
#       check_path(i, verbose = 0)
#       hist_mark = read_file(i)
#       colnames(hist_mark) = tolower(colnames(hist_mark))
#       req_list = c("chr","start","end","state")
#       missing_cols = setdiff(req_list, colnames(hist_mark))
#       if (length(missing_cols) > 0) {
#         stop("Required columns missing in ", i, ": ",paste(missing_cols, collapse = ", "))
#       }
#     }
#     # message("Histone mark instructions columns present.")
#   }
  
#   if (tolower(args[["hic_analysis"]]) == "y") {
#     hic_instruct_dir = args[["hic_instruct_dir"]]
#     check_path(hic_instruct_dir)
#     hic_instruct = read_file(hic_instruct_dir, nrows = -1)
#     colnames(hic_instruct) = tolower(colnames(hic_instruct))
#     req_list = c("hic_dir","genes_present","cell_sorted","atac_cell_types","hic_cell_types")
#     missing_cols = setdiff(req_list, colnames(hic_instruct))
#     if (length(missing_cols) > 0) {
#       stop("Required columns missing in HI-C instructions file: ", paste(missing_cols, collapse = ", "))
#     }
    
#     # check if contents of hic instruct directories exist
#     # for (i in hic_instruct$hic_dir){
#     #   check_path(i, verbose = 0)
#     #   hic_file = read_file(i)
#     #   colnames(hic_file) = tolower(colnames(hic_file))
#     #   req_list <- c("chr1", "start1", "end1", "chr2", "start2", "end2")
#     #   missing_cols = setdiff(req_list, colnames(hic_file))
#     #   if (length(missing_cols) > 0) {
#     #     stop("Required columns missing in ", i, ": ",paste(missing_cols, collapse = ", "))
#     #   }
#     # }
#   }
  
  
#   tf_setting = args[["tf_expr_analysis"]]
#   tf_setting = tolower(tf_setting)
#   tf_list = c("atac","rna","both","none")
#   if (!(tf_setting %in% tf_list)) {
#     stop("Invalid TF option: '", tf_setting, 
#          "'. Allowed values are: ", paste(tf_list, collapse = ", "))
#   }
  
#   if (tolower(args[["tf_expr_analysis"]]) %in% c("rna", "both")) {
#     scrna_instruct_dir = args[["scrna_instruct_dir"]]
#     check_path(scrna_instruct_dir)
#     scrna_instruct = read_file(scrna_instruct_dir, nrows=-1)
#     colnames(scrna_instruct) = tolower(colnames(scrna_instruct))
#     req_list = c("scrna_dir","atac_cell_types","rna_cell_types","one_cell_type_seurat","matrix_loc")
#     missing_cols = setdiff(req_list, colnames(scrna_instruct))
#     if (length(missing_cols) > 0) {
#       stop("Required columns missing in scrna instructions file: ", paste(missing_cols, collapse = ", "))
#     }
    
#     # check if contents of eqtl instruct directories exist
#     for (i in scrna_instruct$scrna_dir){
#       check_path(i, verbose = 0)
#     }
#     # message("scrna instructions columns present.")
#   } 
# }

