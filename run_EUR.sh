#!/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=2         
#SBATCH --mem=20G                 
#SBATCH --time=20:00:00            

# Load necessary modules (if module system is used)
# chmod +x RegSCOUT_main.sh
# module load StdEnv/2023
# module load r/4.4.0

# Call your existing Bash script with all required arguments
./RegSCOUT_main.sh \
  --mode peak_table \
  --finemap Y \
  --SNP_ref /home/ubunkun/Lab/RA_project/RegSCOUT/EUR/input/ \
  --sum_stats /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/GCST90132223_buildGRCh37.tsv_p.tsv \
  --lead_snps /home/ubunkun/Lab/RA_project/RegSCOUT/1preprocess/lead_snps_EUR.txt \
  --plink2_bin /home/ubunkun/anaconda3/envs/bio-R/bin/plink2 \
  --sample_num 97173 \
  --locus_region 1000000 \
  --LD_thr 0.25 \
  --Population EUR \
  --fgwas_src /home/ubunkun/Lab/RA_project/oldpipe/2_GWAS_finemapping/fgwas-0.3.6/src/fgwas \
  --CI_thr 0.95 \
  --output_dir /home/ubunkun/Lab/RA_project/RegSCOUT/EUR/output/ \
  --jaspar_mtx /home/ubunkun/Lab/RA_project/RegSCOUT/common_input/combined_755_motifs_jaspar2024.txt \
  --coaccess_th 0.05 \
  --genome_built hg38 \
  --gencode_dir /home/ubunkun/Lab/RA_project/RegSCOUT/common_input/gencode.v48.annotation.gff3 \
  --peak_th 0.1 \
  --prom_th_up 2000 \
  --prom_th_down 2000 \
  --cic_genomic_window 2000000 \
  --histone_mark_analysis N \
  --hist_mark_instruct_dir "" \
  --tf_expr_analysis atac \
  --scrna_instruct_dir "" \
  --tf_rna_quantile_th 0.25 \
  --hic_eqtl_analysis Y \
  --hic_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/hic_instructions.xlsx \
  --eqtl_instruct_dir /home/ubunkun/Lab/RA_project/RegSCOUT/eqtl_instructions.xlsx \
  --gene_ppa_th 0.05 \


