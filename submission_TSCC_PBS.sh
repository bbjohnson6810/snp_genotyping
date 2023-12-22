#!/bin/bash

#### this is the master submission script for the Palmer lab HS rats genotyping pipeline
#### usage: bash submission_TSCC_PBS.sh

#### set up pipeline arguments and software locations here before executing any code below
#### set the name of the current genotyping round, to be used in downstream outputs
#### set your haplotype dosage preference: Do you want STITCH to impute haplotype dosages?
#### (If so, dosage=yes. If not, dosage=no)

pipeline_arguments=/projects/ps-palmer/hs_rats/round10_2/code_rn7/pipeline_arguments
software=/projects/ps-palmer/hs_rats/round10_2/code_rn7/software/software_location

# parent directory
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
# code (working) directory
code=$(head -n 11 ${pipeline_arguments} | tail -n 1)
# inputs directory
inputs=$(head -n 13 ${pipeline_arguments} | tail -n 1)
# reference genome
ref_gen=$(head -n 15 ${pipeline_arguments} | tail -n 1 | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
# sample metadata produced before the pipeline begins
input_metadata=$(head -n 20 ${pipeline_arguments} | tail -n 1)
# current genotyping round
genotyping_round=$(head -n 27 ${pipeline_arguments} | tail -n 1)

####
# do not edit: this is the formatted metadata file produced during step 1
current_metadata=${dir_path}/demux/sample_sheet.csv


####################################################################################
####################################################################################

echo ""
echo "----------------------     Step 1: Preparation     ----------------------"
chmod u+x ${code}/genotyping/step1_prep.sh
${code}/genotyping/step1_prep.sh ${pipeline_arguments}


####################################################################################
####################################################################################

################# ddGBS sequence demultiplex, trim, alignment #################

# echo ""
# echo "----------------------     Step 2: ddGBS Demultiplex     ----------------------"
# #### a customized command to extract the number of library from metadata
# source activate hs_rats
# num_lib_py=$(cat <<'EOF'
# import pandas as pd
# import sys
# metadata = pd.read_csv(sys.argv[1], dtype=str)
# sys.stdout.write(str(len(metadata["Library_ID"].unique())))
# EOF
# )
# num_lib() { python3 -c "${num_lib_py}" "$@"; }

# #### Overwrite 'current_metadata' with the newly formatted sample_sheet
# current_metadata=${dir_path}/demux/sample_sheet.csv

# #### submit demux array jobs based on the number of libraries, demultiplex per library
# ppn=12
# num_jobs=$(num_lib ${current_metadata})
# STEP2_DEMUX_JOB_ARR=$(qsub -q condo -N demux_kn23 -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
# 						-j n -k oe -m ae \
# 						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
# 						${code}/ddGBS/step2_demux_array_jobs.sh)
# echo "step2_demux: ${STEP2_DEMUX_JOB_ARR}"
# STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )


####################################################################################
####################################################################################

# echo ""
# echo "----------------------     Step 3: ddGBS Alignment     ----------------------"
# #### a customized command to extract the number of sample from metadata
# source activate hs_rats
# num_sample_py=$(cat <<'EOF'
# import pandas as pd
# import sys
# metadata = pd.read_csv(sys.argv[1], dtype=str)
# sys.stdout.write(str(len(metadata["Sample_ID"].unique())))
# EOF
# )
# num_sample() { python3 -c "${num_sample_py}" "$@"; }

# #### submit mapping array jobs based on the number of samples
# ppn=12
# num_jobs=$(num_sample ${current_metadata})
# STEP3_ALIGNMENT_JOB_ARR=$(qsub -q condo -N mapping -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
# 							-j oe -k oe -m ae \
# 							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
# 							-W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
# 							${code}/ddGBS/step3_alignment_array_jobs.sh)
# echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
# STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )


####################################################################################
####################################################################################

################# LcWGS sequence demultiplex, trim, alignment #################

echo ""
echo "----------------------     Step 2: lcWGS Demultiplex     ----------------------"
echo ""
#### a customized command to extract the number of library from metadata
source activate hs_rats
num_lib_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(len(metadata["Library_ID"].unique())))
EOF
)
num_lib() { python3 -c "${num_lib_py}" "$@"; }

## Overwrite 'current_metadata' with the newly formatted sample_sheet
current_metadata=${dir_path}/demux/sample_sheet.csv

## submit demux array jobs based on the number of libraries, demultiplex per library
## -q: the queue to use (condo/hotel)
## -N: name of the job
## -l: job resource list: nodes, cores per node, and walltime to reserve
## -t: integer IDs for the jobs in the array
## -j: output standard output and standard error separately (use -j oe to intermix the two)
## -k: retain both stdout and stderr
## -m: send an email when the job begins, ends, aborts
## -V: export named environment variables into the job environment 
## -v: list of environment variables to export
ppn=12
java_mem=40G
num_jobs=$(num_lib ${current_metadata})
STEP2_DEMUX_JOB_ARR=$(qsub -q condo -N demux_kn23 -l nodes=1:ppn=${ppn},walltime=8:00:00 \
						-t 1-${num_jobs} -j n -k oe -m ea \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",java_mem="${java_mem}" \
						${code}/LcWGS/step2_demux_array_jobs.sh)
echo "step2_demux: ${STEP2_DEMUX_JOB_ARR}"
STEP2_DEMUX_JOB_ARR_id=$(echo "${STEP2_DEMUX_JOB_ARR}" | cut -d '.' -f 1 )


####################################################################################
####################################################################################

echo ""
echo "--------------------     Quality Control 1: MultiQC     --------------------"
echo ""
### submit multiQC array jobs based on the number of libraries

ppn=6
num_jobs=$(num_lib ${current_metadata})
QC1_MULTIQC_JOB=$(qsub -q home -N qc_multiqc_kn23 -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
                       -j n -k oe -m ae  \
                       -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                       -W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
                       ${code}/quality_control/QC1_multiqc_trimming.sh)
echo "QC1_multiQC: ${QC1_MULTIQC_JOB}"


####################################################################################
####################################################################################

echo ""
echo "----------------------     Step 3: lcWGS Alignment     ----------------------"
echo ""
### a customized command to extract the number of sample from metadata
source activate hs_rats
num_sample_py=$(cat <<'EOF'
import pandas as pd
import sys
metadata = pd.read_csv(sys.argv[1], dtype=str)
sys.stdout.write(str(len(metadata["Sample_ID"].unique())))
EOF
)
num_sample() { python3 -c "${num_sample_py}" "$@"; }

### submit mapping array jobs based on the number of samples
ppn=16
java_mem=80G
num_jobs=$(num_sample ${current_metadata})
echo "Submitting job, waiting for resources..."; echo ""
STEP3_ALIGNMENT_JOB_ARR=$(qsub -q condo -N mapping_kn23 -l nodes=1:ppn=${ppn},walltime=8:00:00 -t 1-${num_jobs} \
							-j n -k oe -m ae \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",java_mem="${java_mem}" \
							${code}/LcWGS/step3_alignment_array_jobs.sh)
echo "step3_alignment: ${STEP3_ALIGNMENT_JOB_ARR}"
####							-W depend=afterokarray:${STEP2_DEMUX_JOB_ARR_id} \
STEP3_ALIGNMENT_JOB_ARR_id=$(echo "${STEP3_ALIGNMENT_JOB_ARR}" | cut -d '.' -f 1 )


####################################################################################
####################################################################################

echo ""
echo "--------------------     Quality Control 2: Read Mapping     --------------------"
echo ""
#### submit mapping stats array jobs based on the number of libraries
ppn=16
num_jobs=$(num_sample ${current_metadata})
QC2_MAPPINGRESULT_JOB=$(qsub -q hotel -N qc_mapping -l nodes=1:ppn=${ppn},walltime=24:00:00 -t 1-${num_jobs} \
                             -j n -k oe -m ae  \
                             -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                             ${code}/quality_control/QC2_mappingResult.sh)
####                             -W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
echo "QC2_mapping_stats: ${QC2_MAPPINGRESULT_JOB}"


#####################################################################
######################## Genotyping pipeline ########################
#####################################################################

echo ""
echo "----------------     Step 4a: Genotyping Setup     ----------------"
echo ""

genotyping_round=$(head -n 27 ${pipeline_arguments} | tail -n 1)
prev_metadata=$(head -n 33 ${pipeline_arguments} | tail -n 1)
new_metadata=$(head -n 31 ${pipeline_arguments} | tail -n 1)
genotyping_log=$(head -n 35 ${pipeline_arguments} | tail -n 1)

source activate hs_rats
python3 ${code}/genotyping/util/combine_metadata.py \
    -l ${genotyping_log} \
    -p ${prev_metadata} \
    -m ${new_metadata} \
    -g ${ref_gen} \
    -o ${inputs}/${genotyping_round}
conda deactivate


####################################################################################
####################################################################################

echo ""
echo "----------------     Step 4: Genotype Calling     ----------------"
echo ""

#### change STITCH parameters
k=8 # replace with k for stitch eg. 8
niterations=2 # replace with number of iteration for stitch eg. 2
nGen=100 # replace with nGen for stitch eg. 100
nCore=30 # number of cores for stitch: MUST BE 1 or stitch will crash
method=diploid # replace with method for stitch eg. diploid
tempdir=/oasis/tscc/scratch/bbjohnson # replace with tempdir for stitch eg. /oasis/tscc/scratch/$USER/ - DONE
bamlist=${inputs}/round10_2_bamlist # replace with bamlist for stitch - DONE
sampleNames_file=${inputs}/round10_2_sample_ids # replace with sampleNames for stitch - DONE
pos_dir=/projects/ps-palmer/hs_rats/Ref_panel_mRatBN7.2/STITCH_pos_file # replace with pos_dir for stitch - DONE
reference_panels_loc=$(head -n 29 ${pipeline_arguments} | tail -n 1)

#### submit STITCH variant calling array jobs

ppn=30
source activate hs_rats

declare -A chr_dict
chr_dict=( ['NC_051336.1']=chr1 
			['NC_051337.1']=chr2
			['NC_051338.1']=chr3
			['NC_051339.1']=chr4
			['NC_051340.1']=chr5
			['NC_051341.1']=chr6
			['NC_051342.1']=chr7
			['NC_051343.1']=chr8
			['NC_051344.1']=chr9
			['NC_051345.1']=chr10
			['NC_051346.1']=chr11
			['NC_051347.1']=chr12
			['NC_051348.1']=chr13
			['NC_051349.1']=chr14
			['NC_051350.1']=chr15
			['NC_051351.1']=chr16
			['NC_051352.1']=chr17
			['NC_051353.1']=chr18
			['NC_051354.1']=chr19
			['NC_051355.1']=chr20
			['NC_051356.1']=chrX
			['NC_051357.1']=chrY
			['NC_001665.2']=chrM )

#### loop through chromosomes: one job per chromosome
for chr_nc in "${!chr_dict[@]}"
do	
	chr=${chr_dict[$chr_nc]}
	reference_panels=${reference_panels_loc}_${chr}
	posfile=${pos_dir}/${chr}_STITCH_pos # replace with position file for stitch - DONE
	chunk_file=${dir_path}/${ref_gen}/stitch/${chr}_chunks_${niterations}
	
	#### split the chromosome into chunks
	Rscript ${code}/genotyping/util/STITCH_split_chr.R \
		${posfile}\
		${chunk_file}

	num_chunks=$(cat ${chunk_file} | wc -l)
	((num_jobs=num_chunks-1))

	echo "${chr}: ${num_chunks} chunks"

 	#### run STITCH on each chunk
	STEP4_STITCH_JOB_ARR=$(qsub -q condo -N stitch_${niterations}_${chr}_diploid -l nodes=1:ppn=${ppn}:mem1024,walltime=4:00:00 -t 1-${num_jobs} \
							-j n -k oe -m ae  \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr="${chr_nc}",k="${k}",niterations="${niterations}",nGen="${nGen}",method="${method}",bamlist="${bamlist}",sampleNames_file="${sampleNames_file}",tempdir="${tempdir}",chunk_file="${chunk_file}",nCore=${nCore},posfile="${posfile}",reference_panels="${reference_panels}" \
							${code}/genotyping/step4_stitch_genotypeCalling_array_jobs.sh)
# 	####						-W depend=afterokarray:${STEP3_ALIGNMENT_JOB_ARR_id} \
	echo "step4_${chr}_stitch job array: ${STEP4_STITCH_JOB_ARR}"
	STEP4_STITCH_JOB_ARR_id=$(echo "${STEP4_STITCH_JOB_ARR}" | cut -d '.' -f 1 )
	echo "job ID: ${STEP4_STITCH_JOB_ARR_id}"
	echo ""
done
conda deactivate

# check which STITCH runs successfully completed or failed
${code}/genotyping/util/run_check_stitch_output.sh ${pipeline_arguments}


####################################################################################
####################################################################################

echo ""
echo "----------------     Step 5: Concatenating Variants     ----------------"
echo ""

#### concatenate chunks back together
ppn=30
source activate hs_rats

for chr_nc in "${!chr_dict[@]}"
do	
	chr=${chr_dict[$chr_nc]}
	STEP5_CONCAT_SNPS=$(qsub -q condo -N concat_${chr} -l nodes=1:ppn=${ppn},walltime=4:00:00 \
							-j n -k oe -m ae  \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr="${chr_nc}" \
							${code}/genotyping/step5_concat_variants.sh)
####							-W depend=afterokarray:${STEP4_STITCH_JOB_ARR_id} \
	echo "step5_${chr}_concat: ${STEP5_CONCAT_SNPS}"
	STEP5_CONCAT_SNPS_id=$(echo "${STEP5_CONCAT_SNPS}" | cut -d '.' -f 1 )
done
conda deactivate

####################################################################################
####################################################################################

echo ""
echo "----------------     Step 5b: Rename Chromosomes     ----------------"
echo ""

#### concatenate chunks back together
ppn=6

for chr_nc_val in "${!chr_dict[@]}"
do	
	chr=${chr_dict[$chr_nc_val]}
	chr_nc=${chr_nc_val}
	STEP5_RENAME_CHRS=$(qsub -q condo -N rename_${chr} -l nodes=1:ppn=${ppn},walltime=8:00:00 \
							-j n -k oe -m ae  \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",chr_nc="${chr_nc}",chr="${chr}" \
							${code}/quality_control/new/rename_chrs_each.sh)
####							-W depend=afterokarray:${STEP4_STITCH_JOB_ARR_id} \
	echo "step5_rename_${chr}: ${STEP5_RENAME_CHRS}"
	STEP5_RENAME_CHRS_id=$(echo "${STEP5_RENAME_CHRS}" | cut -d '.' -f 1 )
done

####################################################################################
####################################################################################

echo ""
echo "----------------------     Step 6: Variant Filtering     ----------------------"
echo ""

### file with novel snps to remove from analysis
remove_snps=/projects/ps-palmer/hs_rats/Ref_panel_mRatBN7.2/final_novel_other_exclude_renamed_chr # replace with SNPs position to remove after stitch eg. /projects/ps-palmer/hs_rats/Ref_panel_mRatBN7.2/final_novel_other_exclude
ppn=36
dosage=no

if [ ${dosage} = "yes" ]; then

	STEP6_VARIANT_FILTERING=$(qsub -q hotel -N variant_filtering_HD -l nodes=1:ppn=${ppn},walltime=36:00:00 \
							-j n -k oe -m ae  \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",remove_snps="${remove_snps}" \
							${code}/genotyping/step6_variants_filtering_HD.sh)
	####						-W depend=afterokarray:${STEP5_CONCAT_SNPS_id} \
	echo "step6_variants_filtering_HD: ${STEP6_VARIANT_FILTERING}"
	STEP6_VARIANT_FILTERING_id=$(echo "${STEP6_VARIANT_FILTERING}" | cut -d '.' -f 1 )

elif [ ${dosage} = "no" ]; then

	STEP6_VARIANT_FILTERING=$(qsub -q hotel -N variant_filtering_nonHD -l nodes=1:ppn=${ppn},walltime=30:00:00 \
							-j n -k oe -m ae  \
							-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",remove_snps="${remove_snps}" \
							${code}/genotyping/step6_variants_filtering_nonHD.sh)
	####						-W depend=afterokarray:${STEP5_CONCAT_SNPS_id} \
	echo "step6_variants_filtering_nonHD: ${STEP6_VARIANT_FILTERING}"
	STEP6_VARIANT_FILTERING_id=$(echo "${STEP6_VARIANT_FILTERING}" | cut -d '.' -f 1 )

else echo "Please assign either 'dosage=yes' or 'dosage=no' based on previous STITCH analysis"

fi

####################################################################################
####################################################################################

echo ""
echo "----------------------     Step 7: Mapping & Genotype Summary Statistics     ----------------------"
echo ""

bamlist=${inputs}/${genotyping_round}_bamlist
new_bams=${inputs}/${genotyping_round}_newbams
prev_mapstats=/projects/ps-palmer/hs_rats/round10_1/${ref_gen}/results/mapping_result/mapping_stats_all.csv
genotyping_log=$(head -n 35 ${pipeline_arguments} | tail -n 1)
chr_names=${inputs}/rn7_2_chr_names
chr_nc=${inputs}/rn7_2_chr_nums
raw_stitch_vcf=${dir_path}/genotypes/${genotyping_round}_raw_nofilters.vcf.gz
filtered_stitch_vcf=${dir_path}/${ref_gen}/stitch/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz
n_bams=$(cat ${bamlist} | wc -l)
datestamp=$(date +%y%m%d)
plink_prefix=hs_rats_stitch_n${n_bams}_20${datestamp}

ppn=24

STEP7a_MAPPING_STATS=$(qsub -q hotel -N mapping_stats -l nodes=1:ppn=${ppn},walltime=36:00:00 \
						-j n -k oe -m ae  \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",new_bams="${new_bams}",prev_mapstats="${prev_mapstats}",genotyping_round="${genotyping_round}",chr_names="${chr_names}",chr_nc="${chr_nc}" \
						${code}/quality_control/bbj_step7a_mapping_summary_stats.sh)
echo "step7a_mapping_summary_stats: ${STEP7a_MAPPING_STATS}"
STEP7a_MAPPING_STATS_id=$(echo "${STEP7a_MAPPING_STATS}" | cut -d '.' -f 1 )

STEP7b_GENOTYPE_STATS=$(qsub -q condo -N genotype_stats_info_all -l nodes=1:ppn=${ppn},walltime=8:00:00 \
						-j n -k oe -m ae  \
						-V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}",bamlist="${bamlist}",chr_names="${chr_names}",chr_nc="${chr_nc}",raw_stitch_vcf="${raw_stitch_vcf}",filtered_stitch_vcf="${filtered_stitch_vcf}",plink_prefix="${plink_prefix}",genotyping_round="${genotyping_round}" \
						${code}/quality_control/bbj_step7b_genotype_summary_stats.sh)
echo "step7b_genotype_summary_stats: ${STEP7b_GENOTYPE_STATS}"
STEP7b_GENOTYPE_STATS_id=$(echo "${STEP7b_GENOTYPE_STATS}" | cut -d '.' -f 1 )

####################################################################################
####################################################################################

echo ""
echo "--------------------     Quality Control 3: Genotype Quality     --------------------"
echo ""


### submit genotyping stats array jobs based on the number of library 
ppn=16
new_bams=${dir_path}/${ref_gen}/round10_1_newbams
QC3_GENOTYPERESULT_JOB=$(qsub -q condo -N qc_map_stats -l nodes=1:ppn=${ppn},walltime=8:00:00 \
                              -j n -k oe -m ae  \
                              -V -v pipeline_arguments="${pipeline_arguments}",bamlist="${new_bams}",ppn="${ppn}",software="${software}" \
                              ${code}/quality_control/bbj_QC3.sh)
####                              ${code}/quality_control/QC3_genotypeResult.sh)
####                              -W depend=afterokarray:${STEP6_VARIANT_FILTERING_id} \
echo "QC3_genotype_stats: ${QC3_GENOTYPERESULT_JOB}"



####################################################################################
####################################################################################

echo ""
echo "--------------------     Step 8: Genotype Cleanup & Final Files     --------------------"
echo ""

ppn=24
QC3_GENOTYPERESULT_JOB=$(qsub -q condo -N gtype_cleanup -l nodes=1:ppn=${ppn},walltime=8:00:00 \
                              -j n -k oe -m ae  \
                              -V -v pipeline_arguments="${pipeline_arguments}",ppn="${ppn}",software="${software}" \
                              ${code}/quality_control/bbj_QC3.sh)
####                              ${code}/quality_control/new/gtype_cleanup.sh)
####                              -W depend=afterokarray:${STEP6_VARIANT_FILTERING_id} \
echo "QC3_genotype_stats: ${QC3_GENOTYPERESULT_JOB}"
