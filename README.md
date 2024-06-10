# EffReg_updated
This is the repository for the EffReg pipeline. This pipeline can be used as a stand-alone command line application which integrates the scATAC-seq data with the risk SNPs of diseases to identify the SNP-affected open chromatin sites, Transcription Factors, and Genes. This pipeline can handle both a scATAC-seq data as a Seurat object or processed results of analyzing scATAC-seq data as a cell by peak table along with peak-peak interaction tables. The pipeline has two modes. If a Seurat object is provided, the pipeline can be executed as:
**Bash EffReg_main.sh --mode ATAC_obj --seurat_obj [path to the Seurat object] --Working_dir [Working directory]  --genome_built [genome built of the scATAC-seq data] --genecode_dir  [path to the genecode data of the relevant genome built]  --EffSNP_file [Path to the list of Effect SNPs of the disease]**
And if you want to run the pipeline with cell by peak and peak-peak interaction tables, the command would be:
**Bash EffReg_main.sh --mode peak_table  --Working_dir [Working directory]  --genome_built [genome built of the scATAC-seq data] --genecode_dir  [path to the genecode data of the relevant genome built]  --EffSNP_file [Path to the list of Effect SNPs of the disease]**
In the following, full description of parameters and output files are provided:
## Parameters:
### **Mode:** This option can be used to choose the mode of operation of the pipeline based on input type. If the Seurat object of a scATAC-seq data is provided, this option should be set as ATAC_obj. In this case, the path to the Seurat object should be set with --seurat_obj option. In the case that peak by cell and peak-peak interaction tables are used this option should be set to peak_table. In this case, two sets of tables should be provided in the working directory (set by Working-dir option). The first table should be cell by peak matrix. This table should be exactly named as cell_peak.xlsx and should have three columns named “chr”, “start”, “end” (chromosome, start and end of the peaks), and “cell sub-types” which represents the cell types in which a peak is active. Multiple cells should be comma-separated in this column. **This table should not have any heading and column names should be in the first line.** In addition, cell type-specific peak-peak interaction tables should be provided in the working directory. The file names should be [cell type names].filtered_coaccessible_sites4s.txt. These tables should have three columns exactly named peak1, peak2, and coaccess (coaccessibility value between Peak1 and peak2). 
### **Working_dir:** This option should be used for setting the working directory (the directory address should end with “/”). All the output and intermediate files of the pipeline are generated in this directory. If running the pipeline with processed peak by cell matrix, the required tables of the pipeline should be provided in this directory.
###**Seurat_obj:**If using the pipeline for analyzing Seurat object, this option should be set as the path to the Seurat object.
### **Genome_built:** This option should be used for setting the genome-built of the scATAC-seq data.
###Genecode_dir: This option sets the path to the genecode file which is a reference file with information about genes. Please make sure to use a genecode file compatible with the genome built of your scATAC-seq data. These types of data can be downloaded from the genecode website (https://www.gencodegenes.org).
### **EffSNP_dir:** The path to the Effect SNP table should be provided in this option. Effect SNPs are fine-mapped SNPs of a disease which can change the binding activity of Transcription Factors. This table should have these columns:
SNP: The rs ID of the SNPs should be provided in this column.
CHR: This column represents the chromosome of the SNPs.
Pos: This is the genomic location of the SNPs.
Locus: The locus number of SNPs are provided in this column.
TF: The affected Transcription Factors of the SNPs are shown in this column.
### **Raw_P_value:** This column shows the raw P value of the effect of the SNP on the binding activity of the Transcription Factor.
### **FDR_corr_P_value:** The FDR-corrected P values are provided in this column.

Please make sure that the Effect SNPs have the same genome built as the scATAC-seq data.
