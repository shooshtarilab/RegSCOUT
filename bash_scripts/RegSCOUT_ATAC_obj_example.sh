#!/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=1         
#SBATCH --mem=30G                 
#SBATCH --time=20:00:00            

# Load necessary modules (if module system is used)
chmod +x RegSCOUT_main.sh
module load StdEnv/2023
module load r/4.4.0

# Call your existing Bash script with all required arguments
./RegSCOUT_main.sh \
  --mode ATAC_obj \
  --seurat_obj /project/6043551/rzhan535/scATAC-seq/mouse_brain_E13_5.RData \
  --output_dir /home/rzhan535/scratch/2025/BMI_project/Bella13.5/ \
  --jaspar_mtx /project/6043551/rzhan535/regscout_updated_test/combined_755_motifs_jaspar2024.txt \
  --genome_built hg38 \
  --gencode_dir /project/6043551/rzhan535/gencode.v47.annotation.gff3 \
  --peak_th 0.1 \
  --prom_th_up 2000 \
  --prom_th_down 2000 \