#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
round=$(head -n 27 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
reference_genome=$(head -n 15 ${pipeline_arguments} | tail -n 1)
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch

#### extract software locations from argument files
bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bcftools} = "ERROR" ] || [ ! -f ${bcftools} ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------  HS Rats Genotyping Pipeline  --------------------"
echo "-------------------  Step 6:  Variants filtering  --------------------"
echo "-------------------   (WITH haplotype dosages)    --------------------"
echo "----------------------------------------------------------------------"
echo ""
BEGIN=$(date +%s)

echo "----------------------------------------------------------------------"
echo "------------- Variants filtering with HD using BCFTOOLS  -------------"
echo "----------------------------------------------------------------------"
echo ""
fs=$(ls ${stitch_path}/stitch_*_HD.vcf.gz)

echo "----------------------------------------------------------------------"
echo "Concatenate all stitch HD vcfs:"
echo "${bcftools} concat --threads ${ncpu} --no-version -a -d"
echo "	-O z -o ${stitch_path}/stitch_HD_raw.vcf.gz ${fs}"
echo "----------------------------------------------------------------------"
echo ""

### concatentate all chromosomes into one vcf file (for vcfs with haplotype dosages)
### https://samtools.github.io/bcftools/bcftools.html#concat
###
### NOTE: be sure all files are listed in the correct order - critical for downstream analyses! ####
###
${bcftools} concat --threads ${ncpu} --no-version -a -d none \
	-O z -o ${stitch_path}/stitch_HD_raw.vcf.gz \
    ${stitch_path}/stitch_chr1_HD.vcf.gz \
    ${stitch_path}/stitch_chr2_HD.vcf.gz \
    ${stitch_path}/stitch_chr3_HD.vcf.gz \
    ${stitch_path}/stitch_chr4_HD.vcf.gz \
    ${stitch_path}/stitch_chr5_HD.vcf.gz \
    ${stitch_path}/stitch_chr6_HD.vcf.gz \
    ${stitch_path}/stitch_chr7_HD.vcf.gz \
    ${stitch_path}/stitch_chr8_HD.vcf.gz \
    ${stitch_path}/stitch_chr9_HD.vcf.gz \
    ${stitch_path}/stitch_chr10_HD.vcf.gz \
    ${stitch_path}/stitch_chr11_HD.vcf.gz \
    ${stitch_path}/stitch_chr12_HD.vcf.gz \
    ${stitch_path}/stitch_chr13_HD.vcf.gz \
    ${stitch_path}/stitch_chr14_HD.vcf.gz \
    ${stitch_path}/stitch_chr15_HD.vcf.gz \
    ${stitch_path}/stitch_chr16_HD.vcf.gz \
    ${stitch_path}/stitch_chr17_HD.vcf.gz \
    ${stitch_path}/stitch_chr18_HD.vcf.gz \
    ${stitch_path}/stitch_chr19_HD.vcf.gz \
    ${stitch_path}/stitch_chr20_HD.vcf.gz \
    ${stitch_path}/stitch_chrX_HD.vcf.gz \
    ${stitch_path}/stitch_chrY_HD.vcf.gz \
    ${stitch_path}/stitch_chrM_HD.vcf.gz

END=$(date +%s)
echo "Done concatenating all stitch HD vcfs"
echo "Time elapsed: $(( $END - $BEGIN )) seconds"
echo ""
echo "----------------------------------------------------------------------"
echo "Index the complete HD vcf"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_HD_raw.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""
START=$(date +%s)

#### index the genome-wide vcf file
#### https://samtools.github.io/bcftools/bcftools.html#index
${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_HD_raw.vcf.gz

END=$(date +%s)
echo "Done indexing the complete HD vcf"
echo "Time elapsed: $(( $END - $START )) seconds"
echo ""

filter="INFO_SCORE<0.9"
echo "----------------------------------------------------------------------"
echo "Filter HD vcf by info score"
echo "${bcftools} view --threads ${ppn} -e ${filter}"
echo "	-Oz -o ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz"
echo "	${stitch_path}/stitch_HD_raw.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""
START=$(date +%s)

#### filter out samples with low INFO scores (<0.9)
#### https://samtools.github.io/bcftools/bcftools.html#view
${bcftools} view --threads ${ppn} \
    -e ${filter} \
    -O z -o ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz \
    ${stitch_path}/stitch_HD_raw.vcf.gz

END=$(date +%s)
echo "Done filtering HD vcf by info score"
echo "Time elapsed: $(( $END - $START )) seconds"
echo ""
echo "----------------------------------------------------------------------"
echo "Index the info-filtered vcf"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""
START=$(date +%s)

#### index the filtered vcf file
#### https://samtools.github.io/bcftools/bcftools.html#index
${bcftools} index -t --threads ${ppn} \
    ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz

END=$(date +%s)
echo "Done indexing info-filtered vcf"
echo "Time elapsed: $(( $END - $START )) seconds"
echo ""
echo "----------------------------------------------------------------------"
echo "Filter HD vcf for novel snps"
echo "${bcftools} view --threads ${ppn} -T ^${remove_snps}"
echo "	-O z -o ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz"
echo "	${stitch_path}/stitch_HD_INFO_0.9.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""
START=$(date +%s)

#### remove novel and unknown snps
#### https://samtools.github.io/bcftools/bcftools.html#view
${bcftools} view --threads ${ppn} \
  -T ^${remove_snps} \
  -O z -o ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz \
  ${stitch_path}/stitch_HD_INFO_0.9.vcf.gz

END=$(date +%s)
echo "Done filtering novel snps"
echo "Time elapsed: $(( $END - $START )) seconds"
echo ""
echo "----------------------------------------------------------------------"
echo "Index the novel snp-filtered vcf"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""
START=$(date +%s)

#### index the vcf with removed novel snps
#### https://samtools.github.io/bcftools/bcftools.html#index
${bcftools} index -t --threads ${ppn} \
    ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz

END=$(date +%s)
echo "Done indexing info-filtered vcf"
echo "Time elapsed: $(( $END - $START )) seconds"
echo ""
echo "Concat Variants with HD using BCFTOOLS, total time elapsed: $(( $BEGIN - $START )) seconds"
echo ""

echo "----------------------------------------------------------------------"
echo "----------- Variants filtering with non-HD using BCFTOOLS  -----------"
echo "----------------------------------------------------------------------"
echo ""
START=$(date +%s)

echo "----------------------------------------------------------------------"
echo "Remove HD annotations:"
echo "${bcftools} annotate --threads ${ncpu} -x FORMAT/HD"
echo "	-O z -o ${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""

#### remove haplotype dosage info
#### https://samtools.github.io/bcftools/bcftools.html#annotate
${bcftools} annotate --threads ${ncpu} \
	-x FORMAT/HD \
	-O z -o ${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz ${stitch_path}/stitch_HD_INFO_0.9_rm_discordant.vcf.gz

echo "Done removing HD annotations"
echo ""
echo "----------------------------------------------------------------------"
echo "Index the non-HD vcf"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""

#### index the non-HD vcf
#### https://samtools.github.io/bcftools/bcftools.html#index
${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz

END=$(date +%s)
echo "Done indexing the non-HD vcf"
echo ""
echo "Concat Variants with non-HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"
echo ""



echo "----------------------------------------------------------------------"
echo "------------- Move datasets to final genotype directory)  ------------"
echo "----------------------------------------------------------------------"
echo ""

#### create a new genotypes directory and move final datasets
geno_dir=${dir_path}/genotypes
if [ -d ${geno_dir} ]; then
	echo "folder: ${geno_dir} already exists"
else
	echo "create folder: ${geno_dir}"
	mkdir ${geno_dir}
fi

echo "mv ${stitch_path}/stitch_HD_raw.vcf.gz ${geno_dir}/${round}_HD_raw_nofilters.vcf.gz"

mv ${stitch_path}/stitch_HD_raw.vcf.gz ${geno_dir}/${round}_HD_raw_nofilters.vcf.gz
mv ${stitch_path}/stitch_HD_raw.vcf.gz.tbi ${geno_dir}/${round}_HD_raw_nofilters.vcf.gz.tbi

END=$(date +%s)
echo ""
echo "Done moving vcf files"
echo ""
echo "----------------------------------------------------------------------"
echo "All variant filtering complete!"
echo "Total time elapsed: $(( $BEGIN - $END )) seconds"
echo "----------------------------------------------------------------------"
echo ""


