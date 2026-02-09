#!/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=2         
#SBATCH --mem=60G                 
#SBATCH --time=26:00:00            
#SBATCH --mail-user=jwan2633@uwo.ca
#SBATCH --mail-type=ALL

SCRIPT=~/RegSCOUT

# Load necessary modules (if module system is used)
chmod +x ${SCRIPT}/shooshtari_lab/RegSCOUT_memory.sh
module load StdEnv/2023
module load r/4.4.0

DATA=~/projects/def-pshoosht/jwan2633/regscout_data/inputs

# Call your existing Bash script with all required arguments
${SCRIPT}/bash_scripts/RegSCOUT_main.sh \
  --mode peak_table \
  --finemap Y \
  --snp_ref_dir ${DATA}/ \
  --sum_stats_file ${DATA}/GCST90132223_EUR_buildGRCh37_p.tsv \
  --lead_snps_file ${DATA}/lead_snps_EUR.txt \
  --cicero_dir ${DATA}/cicero/ \
  --gencode_file ${DATA}/gencode.v49lift37.annotation.gff3.gz\
  --peak_cell_file ${DATA}/cell_peak.tsv \
  --sample_num 97173 \
  --population EUR \
  --output_dir ${DATA}/../outputs/ \
  --jaspar_mtx_file ${DATA}/combined_755_motifs_jaspar2024.txt \
  --genome_build hg19 \
  --histone_mark_analysis Y \
  --tf_expr_analysis both \
  --hic_analysis Y \
  --eqtl_analysis Y \
  --instruction_file_dir ${DATA}/ \
  --plink2_dir /cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/plink/2.0.0-a.6.32/bin/plink2 \
  --fgwas_dir /home/jwan2633/projects/def-pshoosht/jwan2633/regscout_data/inputs/fgwas-0.3.6/src/fgwas
~                                                                                                                                                                           
~                                                                                                           