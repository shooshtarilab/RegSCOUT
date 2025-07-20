#!/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=1         
#SBATCH --mem=30G                 
#SBATCH --time=30:00:00            

# Load necessary modules (if module system is used)
chmod +x RegSCOUT_main.sh
module load StdEnv/2023
module load r/4.4.0

# Call your existing Bash script with all required arguments
./RegSCOUT_main.sh \
  --mode ATAC_obj \
  --finemap Y \
  --SNP_ref /project/6043551/rzhan535/regscout_updated_test/snp_ref/ \
  --sum_stats /project/6043551/rzhan535/regscout_updated_test/all_snp_yengo2018_filt.txt \
  --lead_snps /project/6043551/rzhan535/regscout_updated_test/primary_lead_snp_yengo2018.txt \
  --plink2_bin /project/6043551/rzhan535/regscout_updated_test/plink2 \
  --sample_num present \
  --locus_region 1000000 \
  --LD_thr 0.25 \
  --Population EUR \
  --fgwas_src /home/rzhan535/fgwas/src/fgwas \
  --CI_thr 0.95 \
  --seurat_obj /project/6043551/rzhan535/scATAC-seq/mouse_brain_E13_5.RData \
  --output_dir /home/rzhan535/scratch/2025/BMI_project/Bella13.5/ \
  --jaspar_mtx /project/6043551/rzhan535/regscout_updated_test/combined_755_motifs_jaspar2024.txt \
  --coaccess_th 0.05 \
  --genome_built hg38 \
  --gencode_dir /project/6043551/rzhan535/gencode.v47.annotation.gff3 \
  --peak_th 0.1 \
  --prom_th_up 2000 \
  --prom_th_down 2000 \
  --cic_genomic_window 2000000 \
  --histone_mark_analysis Y \
  --hist_mark_instruct_dir /project/6043551/rzhan535/instructions/hist_mark_instructions.xlsx \
  --tf_expr_analysis both \
  --scrna_instruct_dir /project/6043551/rzhan535/instructions/rna_instructions.xlsx \
  --tf_rna_quantile_th 0.25 \
  --hic_eqtl_analysis Y \
  --hic_instruct_dir /project/6043551/rzhan535/instructions/hic_instructions.xlsx \
  --eqtl_instruct_dir /project/6043551/rzhan535/instructions/eqtl_instructions.xlsx \
  --gene_ppa_th 0.05 \

