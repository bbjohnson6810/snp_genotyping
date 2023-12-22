#!/bin/bash

#### read in declared PBS environment variables
pipeline_arguments=$1

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 15 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch


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

for chr_nc in "${!chr_dict[@]}"
do	
	chr=${chr_dict[$chr_nc]}
	pos_chunk_f=${stitch_path}/${chr}_chunks_2
	i=0
	echo "=====================${chr}=====================(ignore the zero)"
	start_pos=0
	while read line; do
		end_pos=${line}
		#### You should re-submit the follow PBSARRAYID for the corresponding chromosome, ignore the Zeros
		if [ ! -f "${stitch_path}/stitch.${chr_nc}.${start_pos}.${end_pos}.vcf.gz" ]; then
		   echo ${i}
		fi
		((i=i+1))
		((start_pos=line+1))
	done < "${pos_chunk_f}"
done
