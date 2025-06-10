#!/bin/bash
set -e
# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --output_dir) output_dir="$2"; shift ;;
    --mode) mode="$2"; shift ;;
    --seurat_obj) seurat_obj="$2"; shift ;;
    --genome_built) genome_built="$2"; shift ;;
    --finemap) finemap="$2"; shift ;;
    --genecode_dir) genecode_dir="$2"; shift;;
    --CI_thr) CI_thr="$2"; shift;;
    --peak_th) peak_th="$2"; shift;;
    --prom_th_down) prom_th_down="$2"; shift;;
    --prom_th_up) prom_th_up="$2"; shift;;
    --SNP_ref) SNP_ref="$2"; shift;;
    --Population) Population="$2"; shift;;
    --sum_stats) sum_stats="$2"; shift;;
    --lead_snps) lead_snps="$2"; shift;;
    --plink2_bin) plink2_bin="$2"; shift;;
    --sample_number) sample_number="$2"; shift;;
    --fgwas_src) fgwas_src="$2"; shift;;
    --ci_gwas_dir) ci_gwas_dir="$2"; shift;;
    
  esac
  shift
done

if [ "$finemap" == "Y" ]; then
    Rscript fgwas_data_prep.R --output_dir "$output_dir" --SNP_ref "$SNP_ref" --Population "$Population" --sum_stats "$sum_stats" --lead_snps "$lead_snps" --plink2_bin "$plink2_bin"
    Rscript fine_map.R --output_dir "$output_dir" --fgwas_src "$fgwas_src" --CI_thr "$CI_thr"
    suff="gwas_CI.txt"
    ci_gwas_dir="${output_dir}${suff}"  
fi
if [ "$mode" == "ATAC_obj" ]; then
    Rscript EffTest_plink.R --afreq_ref "$afreq_ref" --finemapped_gwas "$finemapped_gwas" --CI_thr "$CI_thr" --effect_snp_dir "$output_dir" --genome_built "$genome_built"
    Rscript Peak_cell_extract.R --seurat_obj "$seurat_obj" --output_dir "$output_dir" --peak_th "$peak_th"
    Rscript peak_interaction_extract.R --genome_built "$genome_built" --seurat_obj "$seurat_obj" --output_dir "$output_dir"
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --prom_th_up "$prom_th_up"  --prom_th_down "$prom_th_down"  --genecode_dir "$genecode_dir"

elif [ "$mode" == "peak_table" ]; then
    Rscript EffTest_plink.R --afreq_ref "$afreq_ref" --finemapped_gwas "$finemapped_gwas" --CI_thr "$CI_thr" --effect_snp_dir "$output_dir" --genome_built "$genome_built"
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --prom_th_up "$prom_th_up"  --prom_th_down "$prom_th_down"  --genecode_dir "$genecode_dir"    
fi
