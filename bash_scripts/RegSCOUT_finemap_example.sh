#!/bin/bash

#SBATCH --job-name=RegSCOUT   
#SBATCH --cpus-per-task=1         
#SBATCH --mem=40G                 
#SBATCH --time=14:00:00           

# Load necessary modules (if module system is used)
chmod +x RegSCOUT_main_mod.sh
module load StdEnv/2023
module load r/4.4.0

# Call your existing Bash script with all required arguments
./RegSCOUT_main_mod.sh \
  --finemap Y \
  --output_dir /home/rzhan535/scratch/test_changes_pipeline_2025/ \
  --SNP_ref /project/6043551/rzhan535/regscout_updated_test/snp_ref/ \
  --sum_stats /project/6043551/rzhan535/regscout_updated_test/all_snp_yengo2018_filt.txt \
  --lead_snps /project/6043551/rzhan535/regscout_updated_test/primary_lead_snp_yengo2018.txt \
  --plink2_bin /project/6043551/rzhan535/regscout_updated_test/plink2 \
  --sample_num present \
  --Population EUR \
  --fgwas_src /home/rzhan535/fgwas/src/fgwas \
  --CI_thr 0.95 \
