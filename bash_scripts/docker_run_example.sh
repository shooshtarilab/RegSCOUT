#!/usr/bin/bash


# Enter absolute path of folder containing input files. 
INPUT_DATA=/home/RegSCOUT/RA_IN

# Enter absolute path of local output folder.
OUTPUT_DATA=/home/RegSCOUT/RA_OUT


docker run --rm -u $(id -u):$(id -g) -v "${INPUT_DATA}:${INPUT_DATA}" -v "${OUTPUT_DATA}:${OUTPUT_DATA}" regscout:4.5 \
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
  --instruction_file_dir ${INPUT_DATA}/instruction_files_ubunkun \
