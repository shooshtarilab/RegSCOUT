#!/bin/bash

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --Working_dir) output_dir="$2"; shift ;;
    --mode) mode="$2"; shift ;;
    --seurat_obj) seurat_obj="$2"; shift ;;
    --genome_built) genome_built="$2"; shift ;;
    --EffSNP_file) EffSNP_file="$2"; shift ;;
    --genecode_dir) genecode_dir="$2"; shift;;
  esac
  shift
done

if [ "$mode" == "ATAC_obj" ]; then
    Rscript Peak_cell_extract.R --seurat_obj "$seurat_obj" --output_dir "$output_dir"
    Rscript peak_interaction_extract.R --genome_built "$genome_built" --seurat_obj "$seurat_obj" --output_dir "$output_dir"
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --EffSNP_file "$EffSNP_file" --genecode_dir "$genecode_dir"

elif [ "$mode" == "peak_table" ]; then
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --EffSNP_file "$EffSNP_file" --genecode_dir "$genecode_dir"
    
    
fi
