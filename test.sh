mode="peak_table"
output_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/"
ci_gwas_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/multi_finemap.txt"
genome_built="HG19"
gene_annot_dir="/home/ubunkun/Lab/RA_project/RegSCOUT/inputs/gencode.v48.annotation.gff3.gz"
finemap="Y"
# plink2_bin="/home/ubunkun/anaconda3/envs/bio-R/bin/plnk2"
# fgwas_src="/home/ubunkun/anaconda3/envs/bio-R/bin/fgwas"
sumstats="/home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132224_buildGRCh37_p.tsv"
lead_snp="/home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EAS.txt"
SNP_ref="/home/ubunkun/Lab/RA_project/RegSCOUT/EAS/"
Population="EAS"

Rscript sanity_check.R \
 output_dir=$output_dir \
 ci_gwas_dir=$ci_gwas_dir \
 genome_built=$genome_built \
 gencode_dir=$gene_annot_dir \
 finemap=$finemap \
 sumstats=$sumstats \
 lead_snps=$lead_snps \
 SNP_ref=$SNP_ref
 Population=$Population

# Rscript EffectSNP.R output_dir=$output_dir ci_gwas_dir=$ci_gwas_dir genome_built=$genome_built

#