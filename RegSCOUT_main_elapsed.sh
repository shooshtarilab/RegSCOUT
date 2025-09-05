#!/bin/bash
set -e
set -o pipefail

# Function to run R script and track time
run_rscript() {
    local script_name="$1"
    shift
    
    echo "=== Running $script_name ===" | tee -a "$master_log"
    echo "Start time: $(date)" | tee -a "$master_log"
    
    start_time=$(date +%s)
    Rscript "$script_name" "$@" 2>&1 | tee -a "$master_log"
    end_time=$(date +%s)
    
    elapsed=$((end_time - start_time))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    
    echo "End time: $(date)" | tee -a "$master_log"
    printf "Elapsed time: %02d:%02d:%02d\n" $hours $minutes $seconds | tee -a "$master_log"
    echo "----------------------------------------" | tee -a "$master_log"
}

# Parse named arguments
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --output_dir) output_dir="$2"; shift ;;
    --mode) mode="$2"; shift ;;
    --finemap) finemap="$2"; shift ;;
    --hic_analysis) hic_analysis="$2"; shift;;
    --eqtl_analysis) eqtl_analysis="$2"; shift;;
    --histone_mark_analysis) histone_mark_analysis="$2"; shift;;
    --snp_ref_dir) snp_ref_dir="$2"; shift;;
    --population) population="$2"; shift;;
    --sum_stats_dir) sum_stats_dir="$2"; shift;;
    --lead_snps_dir) lead_snps_dir="$2"; shift;;
    --plink2_dir) plink2_dir="$2"; shift;;
    --sample_num) sample_num="$2"; shift;;
    --locus_region) locus_region="$2"; shift ;;
    --ld_th) ld_th="$2"; shift ;;
    --fgwas_dir) fgwas_dir="$2"; shift;;
    --ci_th) ci_th="$2"; shift;;
    --ci_ppa_th) ci_ppa_th="$2"; shift;;
    --jaspar_mtx_dir) jaspar_mtx_dir="$2"; shift ;;
    --ci_gwas_dir) ci_gwas_dir="$2"; shift;;
    --genome_build) genome_build="$2"; shift ;;
    --seurat_obj_dir) seurat_obj_dir="$2"; shift ;;
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

# Create a logs directory inside output_dir 
ts=$(date +"%Y%m%d_%H%M%S")
master_log="$output_dir/RegSCOUT_pipeline.log"

pipeline_start_time=$(date +%s)
echo "=== Pipeline started at $(date) ===" | tee -a "$master_log"

Rscript sanity_check.R \
    --output_dir "$output_dir" --mode "$mode" --genome_build "$genome_build" \
    --finemap "$finemap" --gencode_dir "$gencode_dir" \
    --snp_ref_dir "$snp_ref_dir" --population "$population" \
    --sum_stats_dir "$sum_stats_dir" --lead_snps_dir "$lead_snps_dir" \
    --plink2_dir "$plink2_dir" --fgwas_dir "$fgwas_dir" \
    --ci_gwas_dir "$ci_gwas_dir" --hic_analysis "$hic_analysis" \
    --eqtl_analysis "$eqtl_analysis" --hic_instruct_dir "$hic_instruct_dir" \
    --eqtl_instruct_dir "$eqtl_instruct_dir" --tf_expr_analysis "$tf_expr_analysis" \
    --scrna_instruct_dir "$scrna_instruct_dir" --hist_mark_instruct_dir "$hist_mark_instruct_dir" \
    --histone_mark_analysis "$histone_mark_analysis" --seurat_obj_dir "$seurat_obj_dir"

if [ "$finemap" == "Y" ]; then
    run_rscript fgwas_data_prep.R \
        --output_dir "$output_dir" --snp_ref_dir "$snp_ref_dir" \
        --population "$population" --sum_stats_dir "$sum_stats_dir" \
        --lead_snps_dir "$lead_snps_dir" --plink2_dir "$plink2_dir" \
        --sample_num "$sample_num" --locus_region "$locus_region" --ld_th "$ld_th"

    run_rscript fine_map.R \
        --output_dir "$output_dir" --fgwas_dir "$fgwas_dir" \
        --ci_th "$ci_th" --ci_ppa_th "$ci_ppa_th"

    suff="gwas_CI.txt"
    ci_gwas_dir="${output_dir}${suff}"
fi

if [ "$mode" == "ATAC_obj" ]; then
    run_rscript EffectSNP.R \
        --output_dir "$output_dir" --jaspar_mtx_dir "$jaspar_mtx_dir" \
        --ci_gwas_dir "$ci_gwas_dir" --genome_build "$genome_build"

    run_rscript Peak_cell_extract.R \
        --seurat_obj_dir "$seurat_obj_dir" --output_dir "$output_dir" \
        --peak_th "$peak_th" --cell_count_th "$cell_count_th"

    run_rscript peak_interaction_extract.R \
        --genome_build "$genome_build" --seurat_obj_dir "$seurat_obj_dir" \
        --output_dir "$output_dir" --coaccess_th "$coaccess_th" \
        --cic_genomic_window "$cic_genomic_window" --peak_th "$peak_th"

    run_rscript Peak_Gene_SNP_Integration.R \
        --output_dir "$output_dir" --prom_th_up "$prom_th_up" \
        --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir"

elif [ "$mode" == "peak_table" ]; then
    run_rscript EffectSNP.R \
        --output_dir "$output_dir" --jaspar_mtx_dir "$jaspar_mtx_dir" \
        --ci_gwas_dir "$ci_gwas_dir" --genome_build "$genome_build"

    run_rscript Peak_Gene_SNP_Integration.R \
        --output_dir "$output_dir" --prom_th_up "$prom_th_up" \
        --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir"
fi

if [ "$tf_expr_analysis" == "atac" ] || [ "$tf_expr_analysis" == "rna" ] || [ "$tf_expr_analysis" == "both" ]; then
    run_rscript TF_expression_analysis.R \
        --output_dir "$output_dir" --prom_th_up "$prom_th_up" \
        --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir" \
        --tf_rna_quantile_th "$tf_rna_quantile_th" --scrna_instruct_dir "$scrna_instruct_dir" \
        --tf_expr_analysis "$tf_expr_analysis"
fi

if [ "$histone_mark_analysis" == "Y" ]; then
    run_rscript histone_mark_analysis.R \
        --output_dir "$output_dir" --hist_mark_instruct_dir "$hist_mark_instruct_dir"
fi

if [ "$hic_analysis" == "Y" ]; then
    run_rscript HiC_Analysis.R \
        --output_dir "$output_dir" --hic_instruct_dir "$hic_instruct_dir" \
        --prom_th_up "$prom_th_up" --prom_th_down "$prom_th_down" --gencode_dir "$gencode_dir"
fi 

if [ "$eqtl_analysis" == "Y" ]; then
    run_rscript eQTL_Analysis.R \
        --output_dir "$output_dir" --eqtl_instruct_dir "$eqtl_instruct_dir" \
        --genome_build "$genome_build"
fi 

if [ -n "$mode" ]; then
    run_rscript final_outputs.R \
        --output_dir "$output_dir" --finemap "$finemap" \
        --ci_gwas_dir "$ci_gwas_dir" --tf_score_th "$tf_score_th" \
        --gene_score_th "$gene_score_th" --gene_sum_ppa_th "$gene_sum_ppa_th"
fi

# Calculate total pipeline time
pipeline_end_time=$(date +%s)
total_elapsed=$((pipeline_end_time - pipeline_start_time))
total_hours=$((total_elapsed / 3600))
total_minutes=$(((total_elapsed % 3600) / 60))
total_seconds=$((total_elapsed % 60))

echo "========================================" | tee -a "$master_log"
echo "=== Pipeline finished at $(date) ===" | tee -a "$master_log"
printf "=== Total pipeline time: %02d:%02d:%02d (hh:mm:ss) ===\n" $total_hours $total_minutes $total_seconds | tee -a "$master_log"
echo "========================================" | tee -a "$master_log"