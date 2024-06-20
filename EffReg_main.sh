#!/bin/bash

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --Working_dir) output_dir="$2"; shift ;;
    --mode) mode="$2"; shift ;;
    --seurat_obj) seurat_obj="$2"; shift ;;
    --genome_built) genome_built="$2"; shift ;;
    --finemapped_gwas) finemapped_gwas="$2"; shift ;;
    --genecode_dir) genecode_dir="$2"; shift;;
    --CI_thr) CI_thr="$2"; shift;;
    --peak_th) peak_th="$2"; shift;;
    --prom_th) prom_th="$2"; shift;;
    --afreq_ref) afreq_ref="$2"; shift;;
    
  esac
  shift
done

if [ "$mode" == "ATAC_obj" ]; then
    rscript EffTest_plink.R --afreq_ref "$afreq_ref" --finemapped_gwas "$finemapped_gwas" --CI_thr "$CI_thr" --effect_snp_dir "$output_dir" --genome_built "$genome_built"
    Rscript Peak_cell_extract.R --seurat_obj "$seurat_obj" --output_dir "$output_dir" --peak_th "$peak_th"
    Rscript peak_interaction_extract.R --genome_built "$genome_built" --seurat_obj "$seurat_obj" --output_dir "$output_dir"
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --prom_th "$prom_th"  --genecode_dir "$genecode_dir"

elif [ "$mode" == "peak_table" ]; then
    rscript EffTest_plink.R --afreq_ref "$afreq_ref" --finemapped_gwas "$finemapped_gwas" --CI_thr "$CI_thr" --effect_snp_dir "$output_di\
r" --genome_built "$genome_built"  
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --prom_th "$prom_th" --genecode_dir "$genecode_dir"
    
    
fi
