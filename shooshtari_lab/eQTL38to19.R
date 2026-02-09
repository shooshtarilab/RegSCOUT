# eQTL hg38 to hg19
# This works for files from eQTL catalogue

suppressPackageStartupMessages(library(liftOver))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(dplyr))
# Use as eQTL38to19.R --output_dir X --eqtl_instruct X --chain_file X
args = commandArgs(trailingOnly = TRUE, asValues = TRUE)

output_dir = args[["output_dir"]]
# output_dir = "/home/ubunkun/Lab/RA_project/RegSCOUT/RA/inputs/eQTL_19/"

eqtl_instruct = args[["eqtl_dir"]]
# eqtl_instruct = "/home/ubunkun/Lab/RA_project/RegSCOUT/RA/instruction_files/eqtl_instructions.tsv"

chain_file = args[["chain_file"]]
# chain_file = "/home/ubunkun/Lab/RA_project/RegSCOUT/shooshtari_lab/hg38ToHg19.over.chain"

chain = import.chain(chain_file)

user_instruct <- read.table(eqtl_instruct, header = T)
num_eqtl = nrow(user_instruct)

 for (i in 1:num_eqtl){
    current_row = user_instruct[i, ]
    eqtl_dir = current_row$eqtl_dir
    tabix = as.logical(current_row$tabix)
    atac_cell_types = unlist(strsplit(current_row$atac_cell_types, split = ","))
    eqtl_d = read.table(file = eqtl_dir, sep = '\t', header = T, comment.char = "")

    gr = makeGRangesFromDataFrame(eqtl_d, 
        keep.extra.columns = TRUE, 
        seqnames.field = "X..chromosome", 
        start.field = "position", 
        end.field = "position")

    new_df = as.data.frame(unlist(liftOver(gr, chain))) %>% dplyr::select(-end, -width, -strand)
    new_df <- new_df[order(new_df$seqnames, new_df$start), ]

    output_file = paste0(output_dir,sub("\\.[[:alnum:]]+$","",basename(current_row$eqtl_dir)))
    write.table(
        new_df,
        file = output_file,
        sep = "\t",
        quote = FALSE,
        row.names = FALSE,
        col.names = FALSE
    )
  
    gz = paste0(output_file,".gz")
    system(paste("bgzip -f", shQuote(output_file)))

    system(paste("tabix -f -s 1 -b 2 -e 2", shQuote(gz)))
 }

