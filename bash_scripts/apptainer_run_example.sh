#!/usr/bin/bash

#SBATCH --job-name=RegSCOUT    
#SBATCH --cpus-per-task=1         
#SBATCH --mem=30G                 
#SBATCH --time=1:00:00            
#SBATCH --output=%x-%j.out

module load StdEnv/2023
module load apptainer/1.4.5

# Enter path of local input folder containing input files.
INPUT_DATA=/home/jwan2633/scratch/RA_IN
# Enter path of local output folder.
OUTPUT_DATA=/home/jwan2633/scratch/RA_OUT
# Enter path of apptainer image.
SIF_IMAGE=/home/jwan2633/scratch/regscout_app/regscout_4.5.2.sif 

# The container will mirr--pwd /RegSCOUTor the path in INPUT_DATA and OUTPUT_DATA inside the container.
apptainer run --pwd /RegSCOUT -B ${INPUT_DATA},${OUTPUT_DATA} ${SIF_IMAGE} \
  --mode peak_table \
  --finemap Y \
  --snp_ref_dir ${INPUT_DATA} \
  --sum_stats_file ${INPUT_DATA}/GCST90132223_EUR_buildGRCh37_p.tsv \
  --lead_snps_file ${INPUT_DATA}/lead_snps_EUR.txt \
  --cicero_dir ${INPUT_DATA}/cicero/ \
  --peak_cell_file ${INPUT_DATA}/cell_peak.tsv \
  --sample_num 97173 \
  --population EUR \
  --output_dir ${OUTPUT_DATA} \
  --jaspar_mtx_file ${INPUT_DATA}/combined_755_motifs_jaspar2024.txt \
  --genome_build hg19 \
  --gencode_file ${INPUT_DATA}/gencode.v49lift37.annotation.gff3.gz  \
  --histone_mark_analysis Y \
  --tf_expr_analysis both \
  --hic_analysis Y \
  --eqtl_analysis Y \
  --instruction_file_dir ${INPUT_DATA}/instruction_files \
