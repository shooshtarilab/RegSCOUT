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
  --finemap Y \
  --snp_ref_dir /home/ubunkun/Lab/RA_project/RegSCOUT/EAS/ \
  --sum_stats_dir /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132224_buildGRCh37_p.tsv \
  --lead_snps_dir /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EAS.txt \
  --sample_num 173633 \
  --population EAS \
  --output_dir /home/ubunkun/Lab/RA_project/RegSCOUT/EAS/ \
  --jaspar_mtx /home/ubunkun/Lab/RA_project/RegSCOUT/inputs/combined_755_motifs_jaspar2024.txt \
  --genome_build hg19 \
  --gencode_dir /home/ubunkun/Lab/RA_project/RegSCOUT/inputs/gencode.v49lift37.annotation.gff3.gz  \
  --histone_mark_analysis N \
  --tf_expr_analysis atac \
  --hic_analysis Y \
  --eqtl_analysis Y \
  --hic_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/hic_instructions.tsv \
  --eqtl_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instructions_spreadsheets/eqtl_instructions.tsv
