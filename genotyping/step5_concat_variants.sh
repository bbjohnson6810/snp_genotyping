#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)

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
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "----------------------  Step 5: Concat Variants  ---------------------"
echo "----------------------------------------------------------------------"
echo ""

# #### USE THIS SECTION ONLY IF ESTIMATING HAPLOTYPE DOSAGES ####
# echo "----------------------------------------------------------------------"
# echo "--------------  Concat Variants with HD using BCFTOOLS  --------------"
# echo "----------------------------------------------------------------------"
# echo ""

# START=$(date +%s)

# fs=$(ls ${stitch_path}/stitch.${chr}.*.vcf.gz)

# echo "----------------------------------------------------------------------"
# echo "HD concat:"
# echo "${bcftools} concat --threads ${ncpu} --no-version -a -d none"
# echo "	-O z -o ${stitch_path}/stitch_${chr}_HD.vcf.gz ${fs}"
# echo "----------------------------------------------------------------------"
# echo ""

# #### save one file per chromosome, including haplotype dosage (HD) info
# #### https://samtools.github.io/bcftools/bcftools.html#concat
# ${bcftools} concat --threads ${ncpu} --no-version -a -d none \
# 	-O z -o ${stitch_path}/stitch_${chr}_HD.vcf.gz ${fs}

# echo "----------------------------------------------------------------------"
# echo "HD index:"
# echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_${chr}_HD.vcf.gz"
# echo "----------------------------------------------------------------------"
# echo ""

# #### index the HD file
# #### https://samtools.github.io/bcftools/bcftools.html#index
# ${bcftools} index -t --threads ${ncpu} \
# 	${stitch_path}/stitch_${chr}_HD.vcf.gz

# END=$(date +%s)
# echo "Concat Variants with HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"


# echo "----------------------------------------------------------------------"
# echo "------------  Concat Variants with non-HD using BCFTOOLS  ------------"
# echo "----------------------------------------------------------------------"
# echo ""
# START=$(date +%s)

# echo "----------------------------------------------------------------------"
# echo "non-HD concat:"
# echo "${bcftools} annotate --threads ${ncpu} -x FORMAT/HD"
# echo "	-O z -o ${stitch_path}/stitch_${chr}_nonHD.vcf.gz ${stitch_path}/stitch_${chr}_HD.vcf.gz"
# echo "----------------------------------------------------------------------"
# echo ""

# #### remove haplotype dosage info and save to a separate file
# #### https://samtools.github.io/bcftools/bcftools.html#annotate
# ${bcftools} annotate --threads ${ncpu} \
# 	-x FORMAT/HD \
# 	-O z -o ${stitch_path}/stitch_${chr}_nonHD.vcf.gz ${stitch_path}/stitch_${chr}_HD.vcf.gz

# echo "----------------------------------------------------------------------"
# echo "non-HD index:"
# echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_${chr}_nonHD.vcf.gz"
# echo "----------------------------------------------------------------------"
# echo ""

# #### index the non-HD file
# #### https://samtools.github.io/bcftools/bcftools.html#index
# ${bcftools} index -t --threads ${ncpu} \
# 	${stitch_path}/stitch_${chr}_nonHD.vcf.gz

# END=$(date +%s)
# echo "Concat Variants with non-HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"
# echo ""
# echo "All concatenation complete!"
# #### END OF SECTION TO USE IF ESTIMATING HAPLOTYPE DOSAGES ####

#### USE THIS SECTION IF YOU ARE **NOT** ESTIMATING HAPLOTYPE DOSAGES ####
echo "----------------------------------------------------------------------"
echo "------------------  Concat Variants using BCFTOOLS  ------------------"
echo "----------------------------------------------------------------------"
echo ""

START=$(date +%s)

fs=$(ls ${stitch_path}/stitch.${chr}.*.vcf.gz)

echo "----------------------------------------------------------------------"
echo "HD concat:"
echo "${bcftools} concat --threads ${ncpu} --no-version -a -d none"
echo "	-O z -o ${stitch_path}/stitch_${chr}_nonHD.vcf.gz ${fs}"
echo "----------------------------------------------------------------------"
echo ""

#### save one file per chromosome, including haplotype dosage (HD) info
#### https://samtools.github.io/bcftools/bcftools.html#concat
${bcftools} concat --threads ${ncpu} --no-version -a -d none \
	-O z -o ${stitch_path}/stitch_${chr}_nonHD.vcf.gz ${fs}

echo "----------------------------------------------------------------------"
echo "HD index:"
echo "${bcftools} index -t --threads ${ncpu} ${stitch_path}/stitch_${chr}_nonHD.vcf.gz"
echo "----------------------------------------------------------------------"
echo ""

#### index the HD file
#### https://samtools.github.io/bcftools/bcftools.html#index
${bcftools} index -t --threads ${ncpu} \
	${stitch_path}/stitch_${chr}_nonHD.vcf.gz

END=$(date +%s)
echo "Concat Variants without HD using BCFTOOLS, time elapsed: $(( $END - $START )) seconds"
