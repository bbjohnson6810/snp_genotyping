#!/bin/bash

#### read in declared PBS environment variables
pipeline_arguments=$1

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
code=$(head -n 11 ${pipeline_arguments} | tail -n 1)
ref_gen=$(head -n 15 ${pipeline_arguments} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
input_metadata=$(head -n 20 ${pipeline_arguments} | tail -n 1)

#### set current directory to be home
cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "--------------------     Step 1: Preparation     ---------------------"
echo "----------------------------------------------------------------------"
echo ""
# echo "----------------------------------------------------------------------"
# echo "--------------------  Construct output directory ---------------------"
# echo "----------------------------------------------------------------------"
# echo ""
# START=$(date +%s)

# #### create the directory structure.
# if [ -d ${dir_path} ]; then
# 	echo "folder: ${dir_path} already exists"
# else
# 	echo "create folder: ${dir_path}"
# 	mkdir ${dir_path}
# fi

# # for new flowcells only, create folders for fastq files
# # if line 24 has no text, create files
# if [ -z $(sed -n "27p" ${pipeline_arguments}) ]; then

#     #### create a folder for demultiplexed fastq files
# 	file=${dir_path}/demux
# 	if [ -d ${file} ]; then
# 		echo "folder: ${file} already exists"
# 	else
# 		echo "create folder: ${file}"
# 		mkdir ${file}
# 		mkdir ${file}/fastq
# 		mkdir ${file}/metrics
# 	fi

# 	#### create a folder for trimmed fastq files
# 	file=${dir_path}/trimmed
# 	if [ -d ${file} ]; then
# 		echo "folder: ${file} already exists"
# 	else
# 		echo "create folder: ${file}"
# 		mkdir ${file}
# 	fi

# 	#### create a folder for sequencing QC results
# 	file=${dir_path}/qc
# 	if [ -d ${file} ]; then
# 		echo "folder: ${file} already exists"
# 	else
# 		echo "create folder: ${file}"
# 		mkdir ${file}
# 	fi

# fi

# #### create a folder to hold input data (metadata, sample lists, bam lists, chr maps)
# file=${dir_path}/inputs
# if [ -d ${file} ]; then
# 	echo "folder: ${file} already exists"
# else
# 	echo "create folder: ${file}"
# 	mkdir ${file}
# fi

# #### make directories to keep qc, sams, bams, stitch, results
# declare -a folders=(${dir_path}/${ref_gen} 
#                     ${dir_path}/${ref_gen}/sams
#                     ${dir_path}/${ref_gen}/bams
#                     ${dir_path}/${ref_gen}/stitch
#                     ${dir_path}/${ref_gen}/results)
# for i in "${folders[@]}"
# do
# 	file=${i}
# 	if [ -d ${file} ]; then
# 		echo "folder: ${file} already exists"
# 		rm -rf ${file}/*
# 	else 
# 		echo "create folder: ${file}"
# 		mkdir ${file}
# 	fi 
# 	#### Creating extra folders for bams
# 	if [ ${file} = ${dir_path}/${ref_gen}/bams ]; then
# 		mkdir ${file}/metrics
# 	fi
# done

# END=$(date +%s)
# echo "Construct output directory, time elapsed: $(( $END - $START )) seconds"


echo "----------------------------------------------------------------------"
echo "----------------- Check and build conda environment ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

#### Check for miniconda
python -m conda
has_conda=$(echo $?)
#### If no miniconda, then download and install it
if [ ${has_conda} != 0 ]; then
	wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
	bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
	rm Miniconda3-latest-Linux-x86_64.sh
fi

#### Check for the 'hs_rats' conda environment
source activate hs_rats
has_conda_env=$(echo $?)
#### If no environment exists, then create it
if [ ${has_conda_env} != 0 ]; then
	~/miniconda3/bin/conda init bash
	source activate base
	# conda create -y -n hs_rats --file ${code}/software/hs_rats_conda_env.yml
	conda env create -y --file ${code}/software/hs_rats_conda_env.yml
	source activate hs_rats
fi

#### Creating from .yml can be finnicky. 
#### If the the environment can't be created due to memory error, try another method
#### NOTE: Be sure this hard-coded version aligns with the desired dependencies in the .yml file! ####
made_conda_env=$(echo $?)
if [ ${made_conda_env} != 0 ]; then
	echo "Retrying 'hs_rats' environment..."
	conda create -y -n hs_rats python=3.10 r-essentials r-stitch r-viridis seaborn multiqc cutadapt -c conda-forge -c bioconda
	source activate hs_rats
fi

echo "Updating hs_rats environment..."
conda update --all -y


END=$(date +%s)
echo "Check and build conda environment, time elapsed: $(( $END - $START )) seconds"


# echo ""
# echo "----------------------------------------------------------------------"
# echo "------------------- Separate metadata by library ---------------------"
# echo "----------------------------------------------------------------------"
# echo ""
# START=$(date +%s)
# #### This part separates and extracts the big sample sheet that Fgbio needs into
# #### several small sample sheets by combination of "pcr_barcode", "library",
# #### and "full_run_id"
# #### This part requires ${input_metadata} to have columns "strain", 
# #### "pcr_barcode", "library_name", "rfid", "project_name", "barcode",
# #### "runid", "fastq_files"
# #### !!!!!!!!!!!!!!!!!!!!!!
# #### The ${code}/genotyping/util/separate_metadata.py probably needs modifications
# #### since original sample sheet (input metadata) always
# #### comes in with DIFFERENT format.
# #### !!!!!!!!!!!!!!!!!!!!!!

# #### This block handles the original sample sheet format.
# if [ ! -f "${dir_path}/demux/sample_sheet.csv" ]; then
# 	#### extract the corresponding sample barcode metadata
# 	source activate hs_rats
# 	python3 ${code}/genotyping/util/separate_metadata.py \
# 		${input_metadata} \
# 		${dir_path}/demux
# 	conda deactivate
# fi

# while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
# 	sleep 60
# done
# END=$(date +%s)
# echo "Separate metadata base on library, time elapsed: $(( $END - $START )) seconds"
