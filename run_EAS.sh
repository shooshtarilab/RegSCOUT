#!/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=1         
#SBATCH --mem=30G                 
#SBATCH --time=20:00:00            

# Load necessary modules (if module system is used)
# chmod +x RegSCOUT_main.sh
# module load StdEnv/2023
# module load r/4.4.0

./RegSCOUT_main.sh \
  --mode peak_table \
  --finemap Y \
  --SNP_ref /home/ubunkun/Lab/RA_project/RegSCOUT/EAS/ \
  --sum_stats /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132224_buildGRCh37_p.tsv \
  --lead_snps /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EAS.txt \
  --plink2_bin /home/ubunkun/anaconda3/envs/bio-R/bin/plink2 \
  --sample_num 173633 \
  --locus_region 1000000 \
  --LD_thr 0.25 \
  --Population EAS \
  --fgwas_src /home/ubunkun/anaconda3/envs/bio-R/bin/fgwas \
  --CI_thr 0.95 \
  --output_dir /home/ubunkun/Lab/RA_project/RegSCOUT/EAS/ \
  --jaspar_mtx /home/ubunkun/Lab/RA_project/RegSCOUT/inputs/combined_755_motifs_jaspar2024.txt\
  --coaccess_th 0.05 \
  --genome_built hg19 \
  --gencode_dir /home/ubunkun/Lab/RA_project/RegSCOUT/inputs/gencode.v48.annotation.gff3.gz \
  --peak_th 0.1 \
  --prom_th_up 2000 \
  --prom_th_down 2000 \
  --cic_genomic_window 2000000 \
  --histone_mark_analysis N \
  --hist_mark_instruct_dir "/home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hist_marks_instructions.xlsx" \
  --tf_expr_analysis atac \
  --scrna_instruct_dir "" \
  --tf_rna_quantile_th 0.25 \
  --hic_eqtl_analysis Y \
  --hic_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hic_instructions.xlsx \
  --eqtl_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/eqtl_instructions.xlsx \
  --gene_ppa_th 0.05 \


# eQTL datasets are in hg38