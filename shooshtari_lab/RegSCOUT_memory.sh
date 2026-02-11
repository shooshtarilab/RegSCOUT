#!/bin/bash
set -e
set -o pipefail

trap 'rm -f "temp_instruct_vars.sh"' EXIT

# Function to run R script and track time
run_rscript() {
    local script_name="$1"
    shift
    
    echo "=== Running $script_name ===" | tee -a "$master_log"
    echo "Start time: $(date)" | tee -a "$master_log"
    
    start_time=$(date +%s)
    /cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/bin/time -v Rscript "$script_name" "$@" 2>&1 | tee -a "$master_log"
    # for Compute Canada replace /usr/bin/time with /cvmfs/soft.computecanada.ca/gentoo/2023/x86-64-v3/usr/bin/time
    end_time=$(date +%s)
    
    elapsed=$((end_time - start_time))
    hours=$((elapsed / 3600))
    minutes=$(((elapsed % 3600) / 60))
    seconds=$((elapsed % 60))
    
    echo "End time: $(date)" | tee -a "$master_log"
    printf "Elapsed time: %02d:%02d:%02d\n" $hours $minutes $seconds | tee -a "$master_log"
    echo "----------------------------------------" | tee -a "$master_log"
}

# Default parameter values
finemap="n" hic_analysis="n" eqtl_analysis="n" histone_mark_analysis="n" tf_expr_analysis="n"
ncores=2

# Parse named arguments
while [[ $# -gt 0 ]]; do
  case $1 in
    --output_dir) output_dir="$2"; shift ;;
    --mode) mode="$2"; shift ;;
    --finemap) finemap="$2"; shift ;;
    --sum_stats_file) sum_stats_file="$2"; shift;;
    --lead_snps_file) lead_snps_file="$2"; shift;;
    --snp_ref_dir) snp_ref_dir="$2"; shift;;
    --population) population="$2"; shift;;
    --sample_num) sample_num="$2"; shift;;
    --locus_region) locus_region="$2"; shift ;;
    --ld_th) ld_th="$2"; shift ;;
    --ci_th) ci_th="$2"; shift;;
    --ci_ppa_th) ci_ppa_th="$2"; shift;;
    --fgwas_dir) fgwas_dir="$2"; shift;;
    --plink2_dir) plink2_dir="$2"; shift;;
    --loci_info_file) loci_info_file="$2"; shift;;
    --ci_gwas_file) ci_gwas_file="$2"; shift;;
    --genome_build) genome_build="$2"; shift ;;
    --gencode_file) gencode_file="$2"; shift;;
    --cicero_dir) cicero_dir="$2"; shift;;
    --peak_cell_file) peak_cell_file="$2"; shift;;
    --jaspar_mtx_file) jaspar_mtx_file="$2"; shift ;;
    --prom_th_down) prom_th_down="$2"; shift;;
    --prom_th_up) prom_th_up="$2"; shift;;
    --seurat_obj_file) seurat_obj_file="$2"; shift ;;
    --cell_count_th) cell_count_th="$2"; shift;;
    --peak_th) peak_th="$2"; shift;;
    --coaccess_th) coaccess_th="$2"; shift;;
    --cic_genomic_window) cic_genomic_window="$2"; shift;;
    --hic_analysis) hic_analysis="$2"; shift;;
    --eqtl_analysis) eqtl_analysis="$2"; shift;;
    --histone_mark_analysis) histone_mark_analysis="$2"; shift;;
    --tf_expr_analysis) tf_expr_analysis="$2"; shift;;
    --instruction_file_dir) instruction_file_dir="$2"; shift;;
    --tf_rna_quantile_th) tf_rna_quantile_th="$2"; shift;;
    --tf_score_th) tf_score_th="$2"; shift;;
    --gene_sum_ppa_th) gene_sum_ppa_th="$2"; shift;;
    --gene_score_th) gene_score_th="$2"; shift;;
    --ncores) ncores="$2"; shift;;
  esac
  shift
done

# Absolute requirements
if [[ -z "${mode:-}" ]]; then
    echo "Error: --mode missing."
    exit 1
fi

if [[ -z "${output_dir:-}" ]]; then
    echo "Error: --output_dir missing."
    exit 1
fi

# Conditional requirements
if [[ "${finemap,,}" == "y" ]]; then
    missing_vars=() 
    for var in snp_ref_dir population sum_stats_file lead_snps_file sample_num; do
        if [[ -z "${!var:-}" ]]; then
            missing_vars+=("--$var")
        fi
    done
    if [[ ${#missing_vars[@]} -gt 0 ]]; then
        echo "Error: The following required parameters are missing for finemapping:"
        echo "  ${missing_vars[*]}"
        exit 1
    fi
fi

if [[ "${mode,,}" == "peak_table" ]]; then
    missing_vars=() 
    for var in genome_build gencode_file peak_cell_file cicero_dir; do
        if [[ -z "${!var:-}" ]]; then
            missing_vars+=("--$var")
        fi
    done
    if [[ ${#missing_vars[@]} -gt 0 ]]; then
        echo "Error: The following required parameters are missing for mode peak_table:"
        echo "  ${missing_vars[*]}"
        exit 1
    fi
fi

if [[ "${finemap,,}" == "n" ]]; then
    if [ -z "$ci_gwas_file" ]; then
        echo "Error: Finemapping data is not provided in --ci_gwas_file"
        exit 1
    fi
fi

if [[ "${mode,,}" == "atac_obj" ]]; then
    if [ -z "$seurat_obj_file" ]; then
        echo "Error: The following required parameter is missing for mode atac_obj:"
        echo "  --seurat_obj_file"
        exit 1
    fi
fi

# Create a logs directory inside output_dir 
ts=$(date +"%Y%m%d_%H%M%S")
master_log="$output_dir/RegSCOUT_pipeline.log"

pipeline_start_time=$(date +%s)
echo "=== Pipeline started at $(date) ===" > "$master_log"

run_rscript sanity_check.R \
    --output_dir "$output_dir" --mode "$mode" --genome_build "$genome_build" \
    --finemap "$finemap" --gencode_file "$gencode_file" \
    --snp_ref_dir "$snp_ref_dir" --population "$population" \
    --sum_stats_file "$sum_stats_file" --lead_snps_file "$lead_snps_file" \
    --plink2_dir "$plink2_dir" --fgwas_dir "$fgwas_dir" --jaspar_mtx_file "$jaspar_mtx_file"\
    --ci_gwas_file "$ci_gwas_file" --peak_cell_file "$peak_cell_file" --cicero_dir "$cicero_dir" --hic_analysis "$hic_analysis" \
    --eqtl_analysis "$eqtl_analysis" --tf_expr_analysis "$tf_expr_analysis" \
    --instruction_file_dir "$instruction_file_dir" --loci_info_file "$loci_info_fil" \
    --histone_mark_analysis "$histone_mark_analysis" --seurat_obj_file "$seurat_obj_file" --ncores "$ncores"

if [[ "${finemap,,}" != "y" && "${loci_info_file,,}" == "" ]]; then
    echo "Warning: Lead and SNP information not set in --loci_info_file. Pipeline will run without this information"
fi

if [ "${finemap,,}" == "y" ]; then
    run_rscript fgwas_data_prep.R \
        --output_dir "$output_dir" --snp_ref_dir "$snp_ref_dir" \
        --population "$population" --sum_stats_file "$sum_stats_file" \
        --lead_snps_file "$lead_snps_file" --plink2_dir "$plink2_dir" \
        --sample_num "$sample_num" --locus_region "$locus_region" --ld_th "$ld_th"

    run_rscript fine_map.R \
        --output_dir "$output_dir" --fgwas_dir "$fgwas_dir" \
        --ci_th "$ci_th" --ci_ppa_th "$ci_ppa_th"

    suff="gwas_CI.txt"
    ci_gwas_file="${output_dir}${suff}"
    loci_info_file="${output_dir}loci_info.txt"
fi

if [ "${mode,,}" == "atac_obj" ]; then
    run_rscript EffectSNP.R \
        --output_dir "$output_dir" --jaspar_mtx_file "$jaspar_mtx_file" \
        --ci_gwas_file "$ci_gwas_file" --genome_build "$genome_build" --ncores "$ncores"

    run_rscript Peak_cell_extract.R \
        --seurat_obj_file "$seurat_obj_file" --output_dir "$output_dir" \
        --peak_th "$peak_th" --cell_count_th "$cell_count_th"
    peak_cell_file="${output_dir}cell_peak.tsv"

    run_rscript peak_interaction_extract.R \
        --genome_build "$genome_build" --seurat_obj_file "$seurat_obj_file" \
        --output_dir "$output_dir" --coaccess_th "$coaccess_th" \
        --cic_genomic_window "$cic_genomic_window" --peak_th "$peak_th"
    cicero_dir="${output_dir}"

    run_rscript Peak_Gene_SNP_Integration.R \
        --output_dir "$output_dir" --prom_th_up "$prom_th_up" \
        --prom_th_down "$prom_th_down" --gencode_file "$gencode_file" \
        --peak_cell_file "$peak_cell_file" --cicero_dir "$cicero_dir"

elif [ "${mode,,}" == "peak_table" ]; then
    run_rscript EffectSNP.R \
        --output_dir "$output_dir" --jaspar_mtx_file "$jaspar_mtx_file" \
        --ci_gwas_file "$ci_gwas_file" --genome_build "$genome_build" --ncores "$ncores"

    run_rscript Peak_Gene_SNP_Integration.R \
        --output_dir "$output_dir" --prom_th_up "$prom_th_up" \
        --prom_th_down "$prom_th_down" --gencode_file "$gencode_file" \
        --peak_cell_file "$peak_cell_file" --cicero_dir "$cicero_dir"
fi

source temp_instruct_vars.sh
rm temp_instruct_vars.sh

if [ "${tf_expr_analysis,,}" == "atac" ] || [ "${tf_expr_analysis,,}" == "rna" ] || [ "${tf_expr_analysis,,}" == "both" ]; then
    run_rscript TF_expression_analysis.R \
        --output_dir "$output_dir" --prom_th_up "$prom_th_up" \
        --prom_th_down "$prom_th_down" \
        --tf_rna_quantile_th "$tf_rna_quantile_th" --scrna_instruct "$scrna_instruct" \
        --tf_expr_analysis "$tf_expr_analysis" --peak_cell_file "$peak_cell_file"
fi

if [ "${histone_mark_analysis,,}" == "y" ]; then
    run_rscript histone_mark_analysis.R \
        --output_dir "$output_dir" --histone_mark_instruct "$histone_mark_instruct"
fi

if [ "${hic_analysis,,}" == "y" ]; then
    run_rscript HiC_Analysis.R \
        --output_dir "$output_dir" --hic_instruct "$hic_instruct" \
        --prom_th_up "$prom_th_up" --prom_th_down "$prom_th_down"
fi 

if [ "${eqtl_analysis,,}" == "y" ]; then
    run_rscript eQTL_Analysis.R \
        --output_dir "$output_dir" --eqtl_instruct "$eqtl_instruct" \
        --genome_build "$genome_build"
fi 

if [ -n "$mode" ]; then
    run_rscript final_outputs.R \
        --output_dir "$output_dir" --finemap "$finemap" \
        --ci_gwas_file "$ci_gwas_file" --tf_score_th "$tf_score_th" \
        --gene_score_th "$gene_score_th" --gene_sum_ppa_th "$gene_sum_ppa_th" \
        --loci_info_file "$loci_info_file"
fi

# Calculate total pipeline time
pipeline_end_time=$(date +%s)
total_elapsed=$((pipeline_end_time - pipeline_start_time))
total_hours=$((total_elapsed / 3600))
total_minutes=$(((total_elapsed % 3600) / 60))
total_seconds=$((total_elapsed % 60))

echo "========================================" | tee -a "$master_log"
echo "=== Pipeline finished at $(date) ===" | tee -a "$master_log"
printf "=== Total pipeline time: %02d:%02d:%02d ===\n" $total_hours $total_minutes $total_seconds | tee -a "$master_log"
echo "========================================" | tee -a "$master_log"