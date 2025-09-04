mode="peak_table"
output_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/"
ci_gwas_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/multi_finemap.txt"
genome_build="HG19"
gencode_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/gencode.v48.annotation.gff3.gz"
finemap="N"
plink2_dir="/home/ubunkun/anaconda3/envs/bio-R/bin/plink2"
fgwas_dir="/home/ubunkun/anaconda3/envs/bio-R/bin/fgwas"
sum_stats_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132224_buildGRCh37_p.tsv"
lead_snps_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EAS.txt"
snp_ref_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/EAS/"
population="EAS"
hic_analysis="Y"
eqtl_analysis="Y"
histone_mark_analysis="Y"
hist_mark_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hist_marks_instructions.tsv"
hic_eqtl_analysis="Y"
hic_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hic_instructions.tsv"
eqtl_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/eqtl_instructions.tsv"
tf_expr_analysis="atac"
scrna_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/scrna_instructions.tsv"

# Rscript sanity_check.R --output_dir "$output_dir" --mode "$mode" --genome_build "$genome_build" --finemap "$finemap" --gencode_dir "$gencode_dir" \
#     --snp_ref_dir "$snp_ref_dir" --population "$population" --sum_stats_dir "$sum_stats_dir" --lead_snps_dir "$lead_snps_dir" --plink2_dir "$plink2_dir" \
#     --fgwas_dir "$fgwas_dir" --ci_gwas_dir "$ci_gwas_dir" --hic_analysis "$hic_analysis" --eqtl_analysis "$eqtl_analysis" --hic_instruct_dir "$hic_instruct_dir" --eqtl_instruct_dir "$eqtl_instruct_dir" \
#     --tf_expr_analysis "$tf_expr_analysis" --scrna_instruct_dir "$scrna_instruct_dir" --hist_mark_instruct_dir "$hist_mark_instruct_dir" \
#     --histone_mark_analysis "$histone_mark_analysis" --seurat_obj_dir "$seurat_obj_dir"


Rscript final_outputs.R --output_dir "$output_dir" --finemap "$finemap" \
--ci_gwas_dir "$ci_gwas_dir" --tf_score_th "$tf_score_th" --gene_score_th \
"$gene_score_th" --gene_sum_ppa_th "$gene_sum_ppa_th"

 


# Rscript EffectSNP.R output_dir=$output_dir ci_gwas_dir=$ci_gwas_dir genome_built=$genome_built

# scrna_instruct_dir hic_instruct_dir eqtl_instruct_dir seurat_obj