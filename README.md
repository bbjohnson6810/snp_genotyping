# snp_genotyping![alt text](https://secureservercdn.net/198.71.233.106/h9j.d46.myftpupload.com/wp-content/uploads/2019/09/palmerlab-logo.png)
![alt text](https://ratgenes.org/wp-content/uploads/2014/11/GWAS_1200x150pxBanner-01.png)

## Source code for genotyping Heterogeneous Stock rats 
### Genotyping Pipeline version 1.0.2 used by the Palmer Lab at the UC San Diego School of Medicine

:information_source: :information_source: :information_source:  **INFOMATION** :information_source: :information_source: :information_source:    
This folder contains all necessary code to impute SNP genotypes. You may want to visit the pipeline flow [here](https://github.com/Palmer-Lab-UCSD/HS-Rats-Genotyping-Pipeline/blob/main/assets/HS_Rats_Lc-WGS_Genotyping_Pipeline_Design.pdf).

## Contents
**[submission_TSCC_PBS.sh](submission_TSCC_PBS.sh)**  
Submission script for the pipeline on  [TSCC](https://www.sdsc.edu/support/user_guides/tscc.html). You may want to change this file accordingly if you are running this pipeline on a different HPC system.  

**[pipeline_arguments](pipeline_arguments)**  
Necessary arguments needed to be set a priori to execute the pipeline

**[genotyping/](genotyping/)**  
All code to process raw Illumina short reads through SNP imputation using STICH

**[quality_control/](quality_control/)**  
All code to check mapping and genotype quality

## Methods

This pipeline imputes biallelic SNPs using both double-digest genotyping-by-sequencing data (ddGBS) and low-coverage whole-genome sequencing data (lcWGS). Briefly, ddGBS sequencing libraries are demultiplexed using [FastX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/commandline.html#fastx_barcode_splitter_usage) and lcWGS libraries are demultiplexed using [fgbio](http://fulcrumgenomics.github.io/fgbio/tools/latest/DemuxFastqs.html). ddGBS sequences are trimmed for quality using [cutadapt](https://cutadapt.readthedocs.io/en/stable/) and lcWGS sequences are trimmed using cutadapt and [bbDuk](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/). All libraries are aligned to the *Rattus norvegicus* [mRatBN7.2](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_015227675.2/) reference genome assembly using [bwa](https://bio-bwa.sourceforge.net/bwa.shtml). SNP genotypes are imputed using [STITCH](https://github.com/rwdavies/STITCH). Imputed genotypes are filtered to keep only those with imputation INFO scores > 0.9. Individual samples missing > 10% of high-INFO genotypes are removed from the dataset, as are samples with heterozygosity rates > 4 or 5 standard deviations from the mean heterozygosity (per library type), and those samples whose phenotypic sex (when available) does not match their genetic sex (determined by read count ratios on the X and Y chromosomes).

## Documentation  
### Before running the pipeline:
Please update the following files to suit your purpose:  
1. Follow the instruction in [software](software) to install required software
2. [pipeline_arguments](pipeline_arguments)
3. [previous_flow_cells_metadata](previous_flow_cells_metadata)
4. [previous_flow_cells_bams](previous_flow_cells_bams)
5. Update the PBS Torque arguments and the corresponding file locations on [submission_TSCC_PBS.sh](submission_TSCC_PBS.sh).  

### Run the pipeline on TSCC:
To run the pipeline on other PBS platform besides TSCC, please change the submission script accordingly.
1. Change the permission of the submission script
```
chmod u+x submission_TSCC_PBS.sh
```
2. Run the submission script
```
./submission_TSCC_PBS.sh
```  
