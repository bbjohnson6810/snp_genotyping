#!/bin/bash

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
code=$(head -n 11 ${pipeline_arguments} | tail -n 1)
inputs=$(head -n 13 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 15 ${pipeline_arguments} | tail -n 1)
ncpu=${ppn}

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch
results_dir=${dir_path}/results
mapping_data=${dir_path}/${ref_gen}/results/mapping_result
mapping_result=${results_dir}/mapping_result
mapping_stats_all=${mapping_data}/${genotyping_round}_mapping_stats_all.csv
# mapping_stats_new=${mapping_data}/${genotyping_round}_mapping_stats_new.csv
mapping_stats_new=${mapping_data}/mapping_stats_new.csv


echo "-------------------------------------------------------------------------"
echo "--------------------- HS Rats Genotyping Pipeline -----------------------"
echo "---------------   Step 7a: Mapping Summary Statistics   -----------------"
echo "-------------------------------------------------------------------------"
echo ""
date
echo ""
BEGIN=$(date +%s)

# #### extract software locations from argument files
# bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# plink2=$(awk 'BEGIN {count = 0} {if ($1 == "Plink2") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# plink1_9=$(awk 'BEGIN {count = 0} {if ($1 == "Plink") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# samtools=$(awk 'BEGIN {count = 0} {if ($1 == "Samtools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
# if [ ${bcftools} = "ERROR" ] || [ ${plink2} = "ERROR" ] || [ ${plink1_9} = "ERROR" ] || [ ${samtools} = "ERROR" ] || [ ! -f ${bcftools} ] || [ ! -f ${plink2} ] || [ ! -f ${plink1_9} ] || [ ! -f ${samtools} ]; then
# 	echo "Error: software_location" 
# 	exit 1
# fi

# #### construct output directories for mapping QC results
# if [ -d ${results_dir} ]; then
# 	echo "QC results folder: ${results_dir} already exists"
# else
# 	echo "Create QC results folder: ${results_dir}"
# 	mkdir ${results_dir}
# fi

# if [ -d ${mapping_data} ]; then
# 	echo "Data folder: ${mapping_data} already exists"
# else
# 	echo "Create intermediate data folder: ${mapping_data}"
# 	mkdir ${mapping_data}
# fi

# if [ -d ${mapping_data}/stats ]; then
#     echo "Mapping stats folder: ${mapping_data}/stats already exists"
# else
#     echo "Create mapping stats folder: ${mapping_data}/stats"
#     mkdir ${mapping_data}/stats
# fi

# if [ -d ${mapping_result} ]; then
# 	echo "Results folder: ${mapping_result} already exists"
# else
# 	echo "Create results folder: ${mapping_result}"
# 	mkdir ${mapping_result}
# fi

# #### create new output files
# if [ -f "${mapping_stats_all}" ]; then
# 	echo "Mapping statistics file: ${mapping_stats_all} already exists"
# 	echo "Note: Output will be appended to the end of the original file"
# else
# 	echo "Create mapping statistics file: ${mapping_stats_all}"
# 	touch ${mapping_stats_all}
# fi

# if [ -f "${mapping_stats_new}" ]; then
# 	echo "Mapping statistics file: ${mapping_stats_new} already exists"
# 	echo "Note: Output will be appended to the end of the original file"
# else
# 	echo "Create mapping statistics file: ${mapping_stats_new}"
# 	touch ${mapping_stats_new}
# fi


# echo "-------------------------------------------------------------------------"
# echo "-------------- Extract mapping statistics for new samples ---------------"
# echo "-------------------------------------------------------------------------"
# echo ""
# START=$(date +%s)

# #### count the number of new bam files associated with the most recent genotyping round
# n_bams=$(cat ${bamlist} | wc -l)

# #### get mapping statistics for each sample
source activate hs_rats
# cnt=0
# for (( line = 1; line <= ${n_bams}; ++line )); do
# 	(( cnt += 1 ))

#   	# save the sample ID
# 	bamfile=$(head -n ${line} ${bamlist} | tail -n 1)
# 	prefix=$(echo ${bamfile} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
# 	sample_id="${prefix%%_mkDup}"
# 	sample_id="${sample_id%%_sorted}"
# 	sample_id="${sample_id%%_trimmed}"
# 	sample_id="${sample_id%%_sorted}"
# 	library_id=$(echo "$sample_id" | sed 's/\(.*\)_[^_]*$/\1/')
# 	rfid=$(echo "$sample_id" | awk -F '_' '{print $NF}')

# 	# count genome-wide mapped reads
# 	# http://www.htslib.org/doc/samtools-flagstat.html
# 	${samtools} flagstat -@ ${ncpu} ${bamfile} > ${mapping_data}/stats/${sample_id}.flagstat &
# 	flagstat=${mapping_data}/stats/${sample_id}.flagstat

# 	# count chromosome-specific mapped reads
# 	# http://www.htslib.org/doc/samtools-idxstats.html
# 	${samtools} idxstats ${bamfile} > ${mapping_data}/stats/${sample_id}.idxstats &
# 	idxstats=${mapping_data}/stats/${sample_id}.idxstats

# 	# get mapping depth statistics
# 	# https://github.com/brentp/mosdepth
# 	# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6030888/
# 	mosdepth --threads ${ncpu} --no-per-base --fasta ${reference_genome} --fast-mode ${mapping_data}/stats/${sample_id} ${bamfile}
# 	mosdepth_summary=${mapping_data}/stats/${sample_id}.mosdepth.summary.txt
	
# 	# organize mapped reads stats
# 	total_reads=$(awk -F" " '{print $1}' ${flagstat} | head -n 1 | tail -n 1)
# 	primary_reads=$(awk -F" " '{print $1}' ${flagstat} | head -n 2 | tail -n 1)
# 	secondary_reads=$(awk -F" " '{print $1}' ${flagstat} | head -n 3 | tail -n 1)
# 	duplicate_primary_reads=$(awk -F" " '{print $1}' ${flagstat} | head -n 6 | tail -n 1)
# 	total_mapped=$(awk -F" " '{print $1}' ${flagstat} | head -n 7 | tail -n 1)
# 	pct_total_mapped=$(awk -F" " '{print $5}' ${flagstat} | head -n 7 | tail -n 1)
# 	pct_total_mapped=${pct_total_mapped#\(} ; pct_total_mapped=${pct_total_mapped%\%}
# 	primary_mapped=$(awk -F" " '{print $1}' ${flagstat} | head -n 8 | tail -n 1)
# 	pct_primary_mapped=$(awk -F" " '{print $6}' ${flagstat} | head -n 8 | tail -n 1)
# 	pct_primary_mapped=${pct_primary_mapped#\(} ; pct_primary_mapped=${pct_primary_mapped%\%}
# 	paired_reads=$(awk -F" " '{print $1}' ${flagstat} | head -n 9 | tail -n 1)
# 	properly_paired=$(awk -F" " '{print $1}' ${flagstat} | head -n 12 | tail -n 1)
# 	pct_properly_paired=$(awk -F" " '{print $6}' ${flagstat} | head -n 12 | tail -n 1)
# 	pct_properly_paired=${pct_properly_paired#\(} ; pct_properly_paired=${pct_properly_paired%\%}
# 	mean_read_depth=$(awk -F "\t" '{if($1=="total") print $4}' ${mosdepth_summary})
# 	min_read_depth=$(awk -F "\t" '{if($1=="total") print $5}' ${mosdepth_summary})
# 	max_read_depth=$(awk -F "\t" '{if($1=="total") print $6}' ${mosdepth_summary})

# 	# organize mapped read stats by chromosome
# 	# saves total mapped reads per indiv per chromosome
# 	for (( chr_line = 1; chr_line <= 23; ++chr_line )); do
# 		chr_name=$(head -n ${chr_line} ${chr_names} | tail -n 1)
# 		chr_num=$(head -n ${chr_line} ${chr_nc} | tail -n 1)
# 		declare "${chr_name}_mapped"=$(awk -v chr_nc=${chr_num} -F "\t" '{if($1==chr_nc) print $3}' ${idxstats})
#         declare "${chr_name}_unmapped"=$(awk -v chr_nc=${chr_num} -F "\t" '{if($1==chr_nc) print $4}' ${idxstats})
#         declare "${chr_name}_mean_depth"=$(awk -v chr_nc=${chr_num} -F "\t" '{if($1==chr_nc) print $4}' ${mosdepth_summary})
#         declare "${chr_name}_min_depth"=$(awk -v chr_nc=${chr_num} -F "\t" '{if($1==chr_nc) print $5}' ${mosdepth_summary})
#         declare "${chr_name}_max_depth"=$(awk -v chr_nc=${chr_num} -F "\t" '{if($1==chr_nc) print $6}' ${mosdepth_summary})
# 	done

 
# 	# append a header to the output file
# 	if [ "${cnt}" == "1" ]; then
		
# 		header_cols=(
# 			"sample_id" "library_id" "rfid" "total_reads" "primary_reads" "secondary_reads" "duplicate_primary_reads" 
# 			"total_mapped" "pct_total_mapped" "primary_mapped" "pct_primary_mapped" "paired_reads" "properly_paired" "mean_read_depth" "min_read_depth" "max_read_depth"
# 			"chr1_mapped" "chr2_mapped" "chr3_mapped" "chr4_mapped" "chr5_mapped" "chr6_mapped" "chr7_mapped" "chr8_mapped"
# 			"chr9_mapped" "chr10_mapped" "chr11_mapped" "chr12_mapped" "chr13_mapped" "chr14_mapped" "chr15_mapped" "chr16_mapped"
# 			"chr17_mapped" "chr18_mapped" "chr19_mapped" "chr20_mapped" "chrX_mapped" "chrY_mapped" "chrM_mapped"
# 			"chr1_unmapped" "chr2_unmapped" "chr3_unmapped" "chr4_unmapped" "chr5_unmapped" "chr6_unmapped" "chr7_unmapped" "chr8_unmapped"
# 			"chr9_unmapped" "chr10_unmapped" "chr11_unmapped" "chr12_unmapped" "chr13_unmapped" "chr14_unmapped" "chr15_unmapped" "chr16_unmapped"
# 			"chr17_unmapped" "chr18_unmapped" "chr19_unmapped" "chr20_unmapped" "chrX_unmapped" "chrY_unmapped" "chrM_unmapped"
# 			"chr1_mean_depth" "chr2_mean_depth" "chr3_mean_depth" "chr4_mean_depth" "chr5_mean_depth" "chr6_mean_depth" "chr7_mean_depth" "chr8_mean_depth"
# 			"chr9_mean_depth" "chr10_mean_depth" "chr11_mean_depth" "chr12_mean_depth" "chr13_mean_depth" "chr14_mean_depth" "chr15_mean_depth" "chr16_mean_depth"
# 			"chr17_mean_depth" "chr18_mean_depth" "chr19_mean_depth" "chr20_mean_depth" "chrX_mean_depth" "chrY_mean_depth" "chrM_mean_depth"
# 			"chr1_min_depth" "chr2_min_depth" "chr3_min_depth" "chr4_min_depth" "chr5_min_depth" "chr6_min_depth" "chr7_min_depth" "chr8_min_depth"
# 			"chr9_min_depth" "chr10_min_depth" "chr11_min_depth" "chr12_min_depth" "chr13_min_depth" "chr14_min_depth" "chr15_min_depth" "chr16_min_depth"
# 			"chr17_min_depth" "chr18_min_depth" "chr19_min_depth" "chr20_min_depth" "chrX_min_depth" "chrY_min_depth" "chrM_min_depth"
# 			"chr1_max_depth" "chr2_max_depth" "chr3_max_depth" "chr4_max_depth" "chr5_max_depth" "chr6_max_depth" "chr7_max_depth" "chr8_max_depth"
# 			"chr9_max_depth" "chr10_max_depth" "chr11_max_depth" "chr12_max_depth" "chr13_max_depth" "chr14_max_depth" "chr15_max_depth" "chr16_max_depth"
# 			"chr17_max_depth" "chr18_max_depth" "chr19_max_depth" "chr20_max_depth" "chrX_max_depth" "chrY_max_depth" "chrM_max_depth"
# 		)

# 		header=$(printf "%s," "${header_cols[@]}")
# 		# remove the trailing comma
# 		header=${header%$','}
# 		# write the header to output files
# 		echo -e "${header}" >> ${mapping_stats_all}
# 		echo -e "${header}" >> ${mapping_stats_new}
# 	fi

# 	# save individual stats to an array
# 	sample_stats=(
# 		${sample_id} ${library_id} ${rfid} ${total_reads} ${primary_reads} ${secondary_reads} ${duplicate_primary_reads} 
# 		${total_mapped} ${pct_total_mapped} ${primary_mapped} ${pct_primary_mapped} ${paired_reads} ${properly_paired} ${mean_read_depth} ${min_read_depth} ${max_read_depth}
# 		${chr1_mapped} ${chr2_mapped} ${chr3_mapped} ${chr4_mapped} ${chr5_mapped} ${chr6_mapped} ${chr7_mapped} ${chr8_mapped}
# 		${chr9_mapped} ${chr10_mapped} ${chr11_mapped} ${chr12_mapped} ${chr13_mapped} ${chr14_mapped} ${chr15_mapped} ${chr16_mapped}
# 		${chr17_mapped} ${chr18_mapped} ${chr19_mapped} ${chr20_mapped} ${chrX_mapped} ${chrY_mapped} ${chrM_mapped}
# 		${chr1_unmapped} ${chr2_unmapped} ${chr3_unmapped} ${chr4_unmapped} ${chr5_unmapped} ${chr6_unmapped} ${chr7_unmapped} ${chr8_unmapped}
# 		${chr9_unmapped} ${chr10_unmapped} ${chr11_unmapped} ${chr12_unmapped} ${chr13_unmapped} ${chr14_unmapped} ${chr15_unmapped} ${chr16_unmapped}
# 		${chr17_unmapped} ${chr18_unmapped} ${chr19_unmapped} ${chr20_unmapped} ${chrX_unmapped} ${chrY_unmapped} ${chrM_unmapped}
# 		${chr1_mean_depth} ${chr2_mean_depth} ${chr3_mean_depth} ${chr4_mean_depth} ${chr5_mean_depth} ${chr6_mean_depth} ${chr7_mean_depth} ${chr8_mean_depth}
# 		${chr9_mean_depth} ${chr10_mean_depth} ${chr11_mean_depth} ${chr12_mean_depth} ${chr13_mean_depth} ${chr14_mean_depth} ${chr15_mean_depth} ${chr16_mean_depth}
# 		${chr17_mean_depth} ${chr18_mean_depth} ${chr19_mean_depth} ${chr20_mean_depth} ${chrX_mean_depth} ${chrY_mean_depth} ${chrM_mean_depth}
# 		${chr1_min_depth} ${chr2_min_depth} ${chr3_min_depth} ${chr4_min_depth} ${chr5_min_depth} ${chr6_min_depth} ${chr7_min_depth} ${chr8_min_depth}
# 		${chr9_min_depth} ${chr10_min_depth} ${chr11_min_depth} ${chr12_min_depth} ${chr13_min_depth} ${chr14_min_depth} ${chr15_min_depth} ${chr16_min_depth}
# 		${chr17_min_depth} ${chr18_min_depth} ${chr19_min_depth} ${chr20_min_depth} ${chrX_min_depth} ${chrY_min_depth} ${chrM_min_depth}
# 		${chr1_max_depth} ${chr2_max_depth} ${chr3_max_depth} ${chr4_max_depth} ${chr5_max_depth} ${chr6_max_depth} ${chr7_max_depth} ${chr8_max_depth}
# 		${chr9_max_depth} ${chr10_max_depth} ${chr11_max_depth} ${chr12_max_depth} ${chr13_max_depth} ${chr14_max_depth} ${chr15_max_depth} ${chr16_max_depth}
# 		${chr17_max_depth} ${chr18_max_depth} ${chr19_max_depth} ${chr20_max_depth} ${chrX_max_depth} ${chrY_max_depth} ${chrM_max_depth}
# 	)

# 	# append stats array to the output file
# 	printf "%s," "${sample_stats[@]}" >> ${mapping_stats_all}
# 	# append a newline character to the file
# 	echo >> ${mapping_stats_all}
# 	# remove the trailing comma
# 	sed -i 's/,$//g' ${mapping_stats_all}

# 	# if a bam file is new to the current sequencing round, also save it to the 'new' output file
# 	if grep -q ${bamfile} ${new_bams}; then
# 		printf "%s," "${sample_stats[@]}" >> ${mapping_stats_new}
# 		# append a newline character to the file
# 		echo >> ${mapping_stats_new}
# 		# remove the trailing comma
# 		sed -i 's/,$//g' ${mapping_stats_new}
# 	fi

# done

echo "mapping_stats_all: ${mapping_stats_all}"
# concatenate old mapstats to new mapstats
python3 ${code}/quality_control/util/concat_mapstats.py \
	-m ${mapping_stats_new} \
	-p ${prev_mapstats} \
	-l ${genotyping_log} \
	-o ${mapping_stats_all}

conda deactivate

FINISH=$(date +%s)
echo "Mapping statistics saved to ${mapping_stats}"
echo ""
echo "-------------------------------------------------------------------------"
echo "-----------------  Mapping summary statistics complete! -----------------"
echo "-------------------------------------------------------------------------"
echo ""
echo "Total time elapsed: $(( $FINISH - $BEGIN )) seconds"
