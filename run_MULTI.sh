#!/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=1         
#SBATCH --mem=30G                 
#SBATCH --time=20:00:00            

# Load necessary modules (if module system is used)
# chmod +x RegSCOUT_main.sh
# module load StdEnv/2023
# module load r/4.4.0

# Call your existing Bash script with all required arguments

  ./RegSCOUT_main.sh \
  --mode peak_table \
  --ci_gwas_dir /home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/multi_finemap.txt \
  --plink2_bin /home/ubunkun/anaconda3/envs/bio-R/bin/plink2 \
  --output_dir  /home/ubunkun/Lab/RA_project/RegSCOUT/MULTI/ \
  --jaspar_mtx /home/ubunkun/Lab/RA_project/RegSCOUT/common_input/combined_755_motifs_jaspar2024.txt \
  --genome_built hg38 \
  --gencode_dir /home/ubunkun/Lab/RA_project/RegSCOUT/common_input/gencode.v48.annotation.gff3 \
  --prom_th_up 2000 \
  --prom_th_down 2000 \
  --histone_mark_analysis N \
  --hist_mark_instruct_dir "" \
  --tf_expr_analysis atac \
  --scrna_instruct_dir "" \
  --tf_rna_quantile_th 0.25 \
  --hic_eqtl_analysis Y \
  --hic_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/hic_instructions_MULTI.xlsx \
  --eqtl_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/eqtl_instructions_MULTI.xlsx \
  --gene_ppa_th 0.05 \

# change EffectSNP cores from 10 to 1

