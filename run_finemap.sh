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
  --output_dir /home/ubunkun/Lab/RA_project/RegSCOUT/EAS/ \
  --snp_ref_dir /home/ubunkun/Lab/RA_project/RegSCOUT/EAS/ \
  --sum_stats_dir /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132224_buildGRCh37_p.tsv \
  --genome_build hg19 \
  --gencode_dir /home/ubunkun/Lab/RA_project/RegSCOUT/inputs/gencode.v19.annotation.gff3.gz \
  --finemap Y \
  --lead_snps_dir /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EAS.txt \
  --sample_num 173633 \
  --population EAS \
  --histone_mark_analysis Y \
  --tf_expr_analysis both \
  --hic_analysis Y \
  --eqtl_analysis Y \
  --hic_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instruction_files/hic_instructions.tsv \
  --eqtl_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instruction_files/eqtl_instructions.tsv \
  --hist_mark_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instruction_files/hist_marks_instructions.tsv \
  --scrna_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/instruction_files/scrna_instructions.tsv \