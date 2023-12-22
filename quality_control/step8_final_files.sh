#!/bin/bash
#SBATCH -J final_files
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -t 24:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node 1
#SBATCH -c 30
#SBATCH --mem-per-cpu 16G
#SBATCH -o sbatch_final_files-%j.o
#SBATCH -e sbatch_final_files-%j.e
#SBATCH --mail-type END,FAIL
#SBATCH --mail-user bbjohnson@health.ucsd.edu
#SBATCH -A csd795

#### arguments included in the qsub command - remove after troubleshooting
pipeline_arguments=/tscc/projects/ps-palmer/hs_rats/round10_2/code_rn7/pipeline_arguments_slurm
bamlist=/tscc/projects/ps-palmer/hs_rats/round10_2/inputs/round10_2_bamlist
new_bams=/tscc/projects/ps-palmer/hs_rats/round10_2/inputs/round10_2_newbams
software=/tscc/projects/ps-palmer/hs_rats/round10_2/code_rn7/software/software_location_slurm
ppn=30

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
code=$(head -n 11 ${pipeline_arguments} | tail -n 1)
inputs=$(head -n 13 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 15 ${pipeline_arguments} | tail -n 1)
geno_round=$(head -n 27 ${pipeline_arguments} | tail -n 1)
n_bams=$(cat ${bamlist} | wc -l)
ncpu=${ppn}

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
util_code=${code}/quality_control/util
metadata=${inputs}/${geno_round}_metadata.csv
genotypes_dir=${dir_path}/genotypes
stitch_path=${dir_path}/${ref_gen}/stitch
stitch_vcf=${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz
results_dir=${dir_path}/results
sample_ids=${inputs}/${geno_round}_sample_ids

#### extract software locations from argument files
bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
plink2=$(awk 'BEGIN {count = 0} {if ($1 == "Plink2") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
plink1_9=$(awk 'BEGIN {count = 0} {if ($1 == "Plink") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
samtools=$(awk 'BEGIN {count = 0} {if ($1 == "Samtools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bcftools} = "ERROR" ] || [ ${plink2} = "ERROR" ] || [ ${plink1_9} = "ERROR" ] || [ ${samtools} = "ERROR" ] || [ ! -f ${bcftools} ] || [ ! -f ${plink2} ] || [ ! -f ${plink1_9} ] || [ ! -f ${samtools} ]; then
	echo "Error: software_location" 
	exit 1
fi


echo "-------------------------------------------------------------------------"
echo "---------------------- HS Rats Genotyping Pipeline ----------------------"
echo "------------    Step 8: Genotype cleanup and final files    -------------"
echo "-------------------------------------------------------------------------"
echo ""
BEGIN=$(date +%s)


echo "------------------------ Genotype cleanup ------------------------"
echo ""

### get the file prefix for the most recently produced genotype dataset
# out_prefix=$(ls ${genotypes_dir}/*.bed  | tail -n 1 | rev |  cut -d '.' -f 2 | cut -d '/' -f 1 | rev)
datestamp=20231213
out_prefix=hs_rats_${geno_round}_n${n_bams}_${datestamp}

## get the file prefix for datasets passing final QC
pass_qc_prefix=$(ls ${results_dir}/${out_prefix}_use_for_analysis*genotype_log.csv | tail -n 1 | cut -d '.' -f 1 | rev | cut -d '/' -f 1 | rev)
pass_qc_prefix="${pass_qc_prefix%%_genotype_log}"



### subset the info-filtered vcf to include only samples that passed QC
### https://samtools.github.io/bcftools/bcftools.html#view
echo "Subset the INFO-filtered vcf for QC-passing samples"
START=$(date +%s)

${bcftools} view -S ${results_dir}/${pass_qc_prefix}_sample_ids --threads ${ncpu} \
	-O z -o ${stitch_path}/${pass_qc_prefix}.vcf.gz \
	${stitch_vcf}

END=$(date +%s)
echo "time elapsed: $(( $END - $START )) seconds"
echo ""

START=$(date +%s)
echo "index ${pass_qc_prefix}.vcf.gz"

${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/${pass_qc_prefix}.vcf.gz

END=$(date +%s)
echo "time elapsed: $(( $END - $START )) seconds"
echo ""

### rename samples in the unfiltered stitch output file
### https://samtools.github.io/bcftools/bcftools.html#reheader
### !! THIS IS THE FINAL GENOTYPE DATASET TO BE USED IN DOWNSTREAM ANALYSES !! ####
START=$(date +%s)
echo "rename samples in ${pass_qc_prefix}.vcf.gz"

${bcftools} reheader -s ${results_dir}/${pass_qc_prefix}_sample_rename \
--threads ${ncpu} -o ${genotypes_dir}/${geno_round}.vcf.gz \
	${stitch_path}/${pass_qc_prefix}.vcf.gz

END=$(date +%s)
echo "time elapsed: $(( $END - $START )) seconds"
echo ""

START=$(date +%s)
echo "index ${geno_round}.vcf.gz"

${bcftools} index -t --threads ${ncpu} \
	${genotypes_dir}/${geno_round}.vcf.gz

END=$(date +%s)
echo "time elapsed: $(( $END - $START )) seconds"
echo ""


echo ""
echo "--------------  Convert STITCH .vcf -> plink binary  ------------------"
echo ""
START=$(date +%s)


### convert STITCH vcf to a plink binary file
### https://www.cog-genomics.org/plink/2.0/input#vcf
### https://www.cog-genomics.org/plink/2.0/input#chr_set
### https://www.cog-genomics.org/plink/2.0/input#allow-extra-chr


# ${plink2} --vcf ${genotypes_dir}/${geno_round}_raw_nofilters.vcf.gz \
# 	--set-missing-var-ids @:# --make-bed --out ${genotypes_dir}/${geno_round}_raw_nofilters

${plink2} --vcf ${genotypes_dir}/${geno_round}.vcf.gz \
	--set-missing-var-ids @:# --make-bed --out ${genotypes_dir}/${geno_round}

END=$(date +%s)
echo "STITCH vcf -> plink time elapsed: $(( $END - $START )) seconds"
echo ""




# echo "------------------------- RMarkdown genotype summary report ------------------------"
# echo ""
# START=$(date +%s)
# source activate hs_rats
# flowcell_ID_py=$(cat <<'EOF'
# import pandas as pd
# import sys
# metadata = pd.read_csv(sys.argv[1], dtype=str)
# sys.stdout.write(str(metadata["runid"].unique()[0]))
# EOF
# )
# flowcell_ID() { python3 -c "${flowcell_ID_py}" "$@"; }

# current_flowcell=$(flowcell_ID ${metadata})

# Rscript ${code}/quality_control/HS_Rats_Genotyping_Summary.r \
# 	${current_flowcell} ${dir_path} ${code} Part1 \
# 	${sex_outliers_Sample_ID}
# conda deactivate 
# END=$(date +%s)
# echo "RMarkdown genotype summary report time elapsed: $(( $END - $START )) seconds"
# echo ""

echo ""
echo "-------------------------------------------------------------------------"
echo "-------------------------------------------------------------------------"
echo ""
echo "-------------------   ALL GENOTYPING IS COMPLETE!   ---------------------"
echo "-------------------------   Enjoy your data!   --------------------------"
echo ""
echo "-------------------------------------------------------------------------"
echo "-------------------------------------------------------------------------"
echo ""
