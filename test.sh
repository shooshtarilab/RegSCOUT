mode="peak_table"
output_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/"
ci_gwas_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/multi_finemap.txt"
genome_built="HG19"
gene_annot_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/gencode.v48.annotation.gff3.gz"
finemap="Y"
# plink2_bin="/home/ubunkun/anaconda3/envs/bio-R/bin/plnk2"
# fgwas_src="/home/ubunkun/anaconda3/envs/bio-R/bin/fgwas"
sum_stats="/home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132224_buildGRCh37_p.tsv"
lead_snp="/home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EAS.txt"
SNP_ref="/home/ubunkun/Lab/RA_project/RegSCOUT/EAS/"
Population="EAS"
hic_eqtl_analysis="Y"
histone_mark_analysis="Y"
hist_mark_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hist_marks_instructions.tsv"
hic_eqtl_analysis="Y"
hic_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hic_instructions.tsv"
eqtl_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/eqtl_instructions.tsv"
tf_expr_analysis="atc"
scrna_instruct_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/scrna_instructions.tsv"

Rscript sanity_check.R \
 output_dir=$output_dir \
 ci_gwas_dir=$ci_gwas_dir \
 genome_built=$genome_built \
 gencode_dir=$gene_annot_dir \
 finemap=$finemap \
 sum_stats=$sum_stats \
 lead_snp=$lead_snp \
 SNP_ref=$SNP_ref \
 Population=$Population \
 histone_mark_analysis=$histone_mark_analysis \
 hist_mark_instruct_dir=$hist_mark_instruct_dir \
 hic_eqtl_analysis=$hic_eqtl_analysis \
 hic_instruct_dir=$hic_instruct_dir \
 eqtl_instruct_dir=$eqtl_instruct_dir \
 tf_expr_analysis=$tf_expr_analysis \
 scrna_instruct_dir=$scrna_instruct_dir
 


 


# Rscript EffectSNP.R output_dir=$output_dir ci_gwas_dir=$ci_gwas_dir genome_built=$genome_built

# scrna_instruct_dir hic_instruct_dir eqtl_instruct_dir seurat_obj