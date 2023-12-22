#!/bin/bash

# cd ${code}
#### extract info from argument files
pipeline_arguments=/tscc/projects/ps-palmer/hs_rats/round10_2/code_rn7/pipeline_arguments_slurm
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
code=$(head -n 11 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 15 ${pipeline_arguments} | tail -n 1)
ncpu=${ppn}
echo "dir_path 1: ${dir_path}"
echo "code 1: ${code}"
echo "reference_genome 1: ${reference_genome}"
echo "ncpu: ${ncpu}"
dir_path=/tscc/projects/ps-palmer/hs_rats/round10_2
code=/tscc/projects/ps-palmer/hs_rats/round10_2/code_rn7
reference_genome=/tscc/projects/ps-palmer/reference_genomes/rat/GCF_015227675.2_mRatBN7.2_genomic.fna
echo "dir_path 2: ${dir_path}"
echo "code 2: ${code}"
echo "reference_genome 2: ${reference_genome}"

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
genotypes_dir=${dir_path}/genotypes
stitch_path=${dir_path}/${ref_gen}/stitch
results_dir=${dir_path}/results
genotype_data=${dir_path}/${ref_gen}/results/genotype_result
genotype_result=${results_dir}/genotype_result
pca_dir=${genotype_result}/pca
echo "ref_gen: ${ref_gen}"
echo "genotypes_dir: ${genotypes_dir}"
echo "stitch_path: ${stitch}"

echo "-------------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ------------------------"
echo "-------------   Step 7b: Genotyping Summary Statistics   ----------------"
echo "-------------------------------------------------------------------------"
echo ""
date
echo ""
BEGIN=$(date +%s)

# ### extract software locations from argument files
# bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# plink2=$(awk 'BEGIN {count = 0} {if ($1 == "Plink2") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# plink1_9=$(awk 'BEGIN {count = 0} {if ($1 == "Plink") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# samtools=$(awk 'BEGIN {count = 0} {if ($1 == "Samtools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# echo "bcftools: ${bcftools}"
# if [ ${bcftools} = "ERROR" ] || [ ${plink2} = "ERROR" ] || [ ${plink1_9} = "ERROR" ] || [ ${samtools} = "ERROR" ] || [ ! -f ${bcftools} ] || [ ! -f ${plink2} ] || [ ! -f ${plink1_9} ] || [ ! -f ${samtools} ]; then
# 	echo "Error: software_location" 
# 	exit 1
# fi
bcftools=/tscc/projects/ps-palmer/software/local/src/bcftools-1.14/bcftools

echo ""
echo "-------------------------------------------------------------------------"
echo "--------------------- Genotyping summary statistics ---------------------"
echo "-------------------------------------------------------------------------"
echo ""

# #### construct output directories for genotyping statistics
# if [ -d ${results_dir} ]; then
# 	echo "QC results folder: ${results_dir} already exists"
# else
# 	echo "Create QC results folder: ${results_dir}"
# 	mkdir ${results_dir}
# fi

# if [ -d ${genotype_result} ]; then
# 	echo "Results folder: ${genotype_result} already exists"
# else
# 	echo "Create results folder: ${genotype_result}"
# 	mkdir ${genotype_result}
# fi

# if [ -d ${pca_dir} ]; then
# 	echo "PCA folder: ${pca_dir} already exists"
# else
# 	echo "Create PCA folder: ${pca_dir}"
# 	mkdir ${pca_dir}
#	mkdir ${pca_dir}/raw
#	mkdir ${pca_dir}/postqc
# fi

# if [ -d ${genotype_data} ]; then
# 	echo "Data folder: ${genotype_data} already exists"
# else
# 	echo "Create intermediate data folder: ${genotype_data}"
# 	mkdir ${genotype_data}
# fi


echo ""
echo "------------------------ Extract STITCH INFO scores ----------------------"
echo ""
START=$(date +%s)

echo "dir_path: ${dir_path}"
echo "code: ${code}"
echo "genotype_data: ${genotype_data}"
echo "ppn: ${ppn}"
echo "pipeline_arguments: ${pipeline_arguments}"
#### create files for snp INFO scores
if [ -f "${genotype_data}/${plink_prefix}_INFO_all_snps" ]; then
	echo "File ${genotype_data}/${plink_prefix}_INFO_all_snps already exists"
else
	echo "Saving INFO scores for all genotyped snps: ${genotype_data}/${plink_prefix}_INFO_all_snps"
	touch ${genotype_data}/${plink_prefix}_INFO_all_snps
	# add a header to the new file
	echo -e "CHROM\tPOS\tINFO_SCORE" >> ${genotype_data}/${plink_prefix}_INFO_all_snps
	# extract info scores
	# https://samtools.github.io/bcftools/bcftools.html#query
	${bcftools} query -f '%CHROM\t%POS\t%INFO/INFO_SCORE\n' \
		${raw_stitch_vcf} >> ${genotype_data}/${plink_prefix}_INFO_all_snps

fi
echo ""

# subset all SNPs with INFO >= 0.9 to a new file 
if [ -f "${genotype_data}/${genotyping_round}_INFO_0.9_snps" ]; then
	echo "File ${genotype_data}/${genotyping_round}_INFO_0.9_snps already exists"
else
	echo "Saving INFO scores for all filtered snps: ${genotype_data}/${genotyping_round}_INFO_0.9_snps"
	awk -F'\t' '$3 > 0.9' ${genotype_data}/${genotyping_round}_INFO_all_snps > ${genotype_data}/${genotyping_round}_INFO_0.9_snps

fi
echo ""

END=$(date +%s)
echo "INFO scores for all SNPs saved in ${genotype_data}/${genotyping_round}_INFO_all_snps"
echo "INFO scores for all filtered SNPs saved in ${genotype_data}/${genotyping_round}_INFO_0.9_snps" 
echo "INFO scores time elapsed: $(( $END - $START )) seconds"
echo ""


# echo ""
# echo "--------------  Convert STITCH .vcf -> plink binary  ------------------"
# echo ""
# START=$(date +%s)


# #### convert STITCH vcf to a plink binary file
# #### https://www.cog-genomics.org/plink/2.0/input#vcf
# #### https://www.cog-genomics.org/plink/2.0/input#chr_set
# #### https://www.cog-genomics.org/plink/2.0/input#allow-extra-chr
# ${plink2} --vcf ${filtered_stitch_vcf} \
# 	--set-missing-var-ids @:# --make-bed --chr-set 20 no-xy  --out ${stitch_path}/${plink_prefix}

# END=$(date +%s)
# echo "STITCH vcf -> plink time elapsed: $(( $END - $START )) seconds"
# echo ""


# echo ""
# echo "----------------  Calculate per-snp missing rates  --------------------"
# echo ""
# START=$(date +%s)

# #### calculate per-snp missingness rates (% of samples with missing data at each locus)
# #### https://www.cog-genomics.org/plink/2.0/basic_stats#missing
# START=$(date +%s)
# ${plink2} --bfile ${stitch_path}/${plink_prefix} \
# 	--missing --chr-set 20 no-xy --out ${genotype_data}/${plink_prefix}

# END=$(date +%s)
# echo "SNP missing rate time elapsed: $(( $END - $START )) seconds"
# echo ""


# echo ""
# echo "------------  Calculate per-sample heterozygosity rates  ------------------"
# echo ""
# START=$(date +%s)

# #### https://www.cog-genomics.org/plink/2.0/basic_stats#het
# #### https://www.cog-genomics.org/plink/2.0/formats#het
# ${plink2} --bfile ${stitch_path}/${plink_prefix} \
# 	--het --chr-set 20 no-xy --out ${genotype_data}/${plink_prefix}

# END=$(date +%s)
# echo "Sample heterozygosity time elapsed: $(( $END - $START )) seconds"
# echo ""


# echo ""
# echo "----------------  Calculate per-snp allele frequencies  ------------------"
# echo ""
# START=$(date +%s)

# #### calculate MAF 
# #### https://www.cog-genomics.org/plink/2.0/basic_stats#freq
# ${plink2} --bfile ${stitch_path}/${plink_prefix} \
# 	--freq --chr-set 20 no-xy --out ${genotype_data}/${plink_prefix}

# END=$(date +%s)
# echo "Allele frequencies time elapsed: $(( $END - $START )) seconds"
# echo ""


# echo ""
# echo "-----------  Calculate per-snp Hardy-Weinberg test statistics  ------------"
# echo ""
# START=$(date +%s)

# #### https://www.cog-genomics.org/plink/2.0/basic_stats#hardy
# ${plink2} --bfile ${stitch_path}/${plink_prefix} \
# 	--hardy --chr-set 20 no-xy --out ${genotype_data}/${plink_prefix}

# END=$(date +%s)
# echo "HWE statistics time elapsed: $(( $END - $START )) seconds"


# echo ""
# echo "----------------------  PCA on all genotypes  --------------------------"
# echo ""
# START=$(date +%s)

### TO DO ###
# just create a GRM for all gtyped samples
# then plot and re-plot by directly subsetting the GRM
# ex: plink --bfile your_data --make-grm --out your_output_prefix
# conduct pca in R: pca_result <- prcomp(cov_mat)
#############

# #### https://www.cog-genomics.org/plink/2.0/strat#pca
# ${plink2} --bfile ${stitch_path}/${plink_prefix} \
# 	--pca --chr-set 20 no-xy -out ${genotype_data}/${plink_prefix}

# END=$(date +%s)
# echo "SNPs PCA time elapsed: $(( $END - $START )) seconds"
# echo ""

FINISH=$(date +%s)
echo ""
echo "-------------------------------------------------------------------------"
echo "----------------  Genotype summary statistics complete! -----------------"
echo "-------------------------------------------------------------------------"
echo ""
echo "Total time elapsed: $(( $FINISH - $BEGIN )) seconds"
