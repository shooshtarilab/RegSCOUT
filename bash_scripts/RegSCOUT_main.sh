#!/bin/bash
set -e
# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --output_dir) output_dir="$2"; shift ;;
    --mode) mode="$2"; shift ;;
    --finemap) finemap="$2"; shift ;;
    --hic_eqtl_analysis) hic_eqtl_analysis="$2"; shift;;
    --histone_mark_analysis) histone_mark_analysis="$2"; shift;;
    --snp_ref) snp_ref="$2"; shift;;
    --population) population="$2"; shift;;
    --sum_stats) sum_stats="$2"; shift;;
    --lead_snps) lead_snps="$2"; shift;;
    --plink2_bin) plink2_bin="$2"; shift;;
    --sample_num) sample_num="$2"; shift;;
    --locus_region) locus_region="$2"; shift ;;
    --ld_th) ld_th="$2"; shift ;;
    --fgwas_src) fgwas_src="$2"; shift;;
    --ci_th) ci_th="$2"; shift;;
    --ci_ppa_th) ci_ppa_th="$2"; shift;;
    --jaspar_mtx) jaspar_mtx="$2"; shift ;;
    --ci_gwas_dir) ci_gwas_dir="$2"; shift;;
    --genome_build) genome_build="$2"; shift ;;
    --seurat_obj) seurat_obj="$2"; shift ;;
    --cell_count_th) cell_count_th="$2"; shift;;
    --peak_th) peak_th="$2"; shift;;
    --coaccess_th) coaccess_th="$2"; shift;;
    --cic_genomic_window) cic_genomic_window="$2"; shift;;
    --prom_th_down) prom_th_down="$2"; shift;;
    --prom_th_up) prom_th_up="$2"; shift;;
    --gencode_dir) gencode_dir="$2"; shift;;
    --hic_instruct_dir) hic_instruct_dir="$2"; shift;;
    --eqtl_instruct_dir) eqtl_instruct_dir="$2"; shift;;
    --tf_expr_analysis) tf_expr_analysis="$2"; shift;;
    --scrna_instruct_dir) scrna_instruct_dir="$2"; shift;;
    --tf_rna_quantile_th) tf_rna_quantile_th="$2"; shift;;
    --hist_mark_instruct_dir) hist_mark_instruct_dir="$2"; shift;;
    --tf_score_th) tf_score_th="$2"; shift;;
    --gene_sum_ppa_th) gene_sum_ppa_th="$2"; shift;;
    --gene_score_th) gene_score_th="$2"; shift;;
  esac
  shift
done

if [ "$finemap" == "Y" ]; then
    Rscript fgwas_data_prep.R --output_dir "$output_dir" --snp_ref "$snp_ref" --population "$population" --sum_stats "$sum_stats" --lead_snps "$lead_snps" --plink2_bin "$plink2_bin" --sample_num "$sample_num" --locus_region "$locus_region" --ld_th "$ld_th"
    Rscript fine_map.R --output_dir "$output_dir" --fgwas_src "$fgwas_src" --ci_th "$ci_th" --ci_ppa_th "$ci_ppa_th"
    suff="gwas_CI.txt"
    ci_gwas_dir="${output_dir}${suff}"  
fi

if [ "$mode" == "ATAC_obj" ]; then
    Rscript EffectSNP.R --output_dir "$output_dir" --jaspar_mtx "$jaspar_mtx" --ci_gwas_dir "$ci_gwas_dir" --genome_build "$genome_build"
    Rscript Peak_cell_extract.R --seurat_obj "$seurat_obj" --output_dir "$output_dir" --peak_th "$peak_th" --cell_count_th "$cell_count_th"
    Rscript peak_interaction_extract.R --genome_build "$genome_build" --seurat_obj "$seurat_obj" --output_dir "$output_dir" --coaccess_th "$coaccess_th" --cic_genomic_window "$cic_genomic_window" --peak_th "$peak_th"
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --prom_th_up "$prom_th_up"  --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir"

elif [ "$mode" == "peak_table" ]; then
    Rscript EffectSNP.R --output_dir "$output_dir" --jaspar_mtx "$jaspar_mtx" --ci_gwas_dir "$ci_gwas_dir" --genome_build "$genome_build"
    Rscript Peak_Gene_SNP_Integration.R --output_dir "$output_dir" --prom_th_up "$prom_th_up"  --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir"    
fi

if [ "$tf_expr_analysis" == "atac" ] || [ "$tf_expr_analysis" == "rna" ] || [ "$tf_expr_analysis" == "both" ]; then
    Rscript TF_expression_analysis.R --output_dir "$output_dir" --prom_th_up "$prom_th_up" --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir" --tf_rna_quantile_th "$tf_rna_quantile_th" --scrna_instruct_dir "$scrna_instruct_dir" --tf_expr_analysis "$tf_expr_analysis"
fi

if [ "$histone_mark_analysis" == "Y" ]; then
    Rscript histone_mark_analysis.R --output_dir "$output_dir" --hist_mark_instruct_dir "$hist_mark_instruct_dir" 
fi

if [ "$hic_eqtl_analysis" == "Y" ]; then
    Rscript HiC_Analysis.R --output_dir "$output_dir" --hic_instruct_dir "$hic_instruct_dir" --prom_th_up "$prom_th_up" --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir"
    Rscript eQTL_Analysis.R --output_dir "$output_dir" --eqtl_instruct_dir "$eqtl_instruct_dir"
    Rscript final_outputs.R --output_dir "$output_dir" --finemap "$finemap" --ci_gwas_dir "$ci_gwas_dir" --tf_score_th "$tf_score_th" --gene_score_th "$gene_score_th" --gene_sum_ppa_th "$gene_sum_ppa_th"

elif [[ -z "$hic_eqtl_analysis" && -z "$mode" ]]; then
    Rscript final_outputs.R --output_dir "$output_dir" --finemap "$finemap" --ci_gwas_dir "$ci_gwas_dir" --tf_score_th "$tf_score_th" --gene_score_th "$gene_score_th" --gene_sum_ppa_th "$gene_sum_ppa_th"
fi
