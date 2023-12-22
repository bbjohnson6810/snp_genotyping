#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
fastq_dir=$(head -n 17 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
demux_dir=${dir_path}/demux
trimmed_dir=${dir_path}/trimmed

#### extract software locations from argument files
java=$(awk 'BEGIN {count = 0} {if ($1 == "Java") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
fgbio=$(awk 'BEGIN {count = 0} {if ($1 == "Fgbio") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
BBMap=$(awk 'BEGIN {count = 0} {if ($1 == "BBMap") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${java} = "ERROR" ] || [ ${fgbio} = "ERROR" ] || [ ${BBMap} = "ERROR" ] || [ ! -f ${java} ] || [ ! -f ${fgbio} ] || [ ! -d ${BBMap} ]; then
	echo "Error: software_location" 
	exit 1
fi

cd $HOME

# echo "----------------------------------------------------------------------"
# echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
# echo "-----------------     Step 2: lcWGS Demultiplex     ------------------"
# echo "----------------------------------------------------------------------"
# echo ""
# echo "----------------------------------------------------------------------"
# echo "----------------------  Demultiplex with Fgbio -----------------------"
# echo "----------------------------------------------------------------------"
# echo ""
# START=$(date +%s)

#### !!!!!!!!!!!!!!!!!!!!!!
#### The following part may need modifications
#### check separate_metadata.py for format.
#### !!!!!!!!!!!!!!!!!!!!!!
sample_sheet=$(ls ${demux_dir}/SampleSheet_*.csv | head -n ${PBS_ARRAYID} | tail -n 1)
flow_cell=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1 | cut -d '_' -f4-)
pre_demux_fastqs=$(head -n 2 ${sample_sheet} | tail -n 1 | cut -d ',' -f6)
pre_demux_fastq_R1=$(cut -d ';' -f1 <<< ${pre_demux_fastqs})
pre_demux_fastq_R2=$(cut -d ';' -f2 <<< ${pre_demux_fastqs} | sed 's/^ *//g')
metrics_name=$(echo ${sample_sheet} | rev | cut -d '/' -f1 | rev | cut -d '.' -f1)


echo "----------------------------------------------------------------------"
echo "${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap "
echo "-jar ${fgbio} DemuxFastqs "
echo "--inputs ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} "
echo "         ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} "
echo "--metadata ${sample_sheet} "
echo "--read-structures 8B12M+T 8M+T "
echo "--output-type=Fastq "
echo "--threads ${ncpu} "
echo "--output ${demux_dir}/fastq "
echo "--metrics ${demux_dir}/metrics/${metrics_name}_demux_barcode_metrics.txt"
echo "----------------------------------------------------------------------"
echo ""

#### check for necessary files 
if [ ! -f "${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1}" ] || [ ! -f "${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2}" ]; then 
	echo "Error: ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} or ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} doesn't exist" 
	exit 1
fi

if [ ! -f "${sample_sheet}" ]; then 
	echo "Error: ${sample_sheet} doesn't exist. Check step1_prep.sh and separate_metadata.py output" 
	exit 1
fi

#### http://fulcrumgenomics.github.io/fgbio/tools/latest/DemuxFastqs.html
#### https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
#### https://fulcrumgenomics.github.io/fgbio/validate-read-structure.html
#### read-structures : two structures used (for 2 paired-end reads)
#### 8B12M+T : identifies an 8bp sample barcode, a 12bp molecular barcode, 
####           remaining bps as template (genomic DNA) from R1
#### 8M+T : identifies an 8bp molecular barcode, remaining bps as template (genomic DNA) from R2


${java} -Xmx${java_mem} -XX:+AggressiveOpts -XX:+AggressiveHeap \
	-jar ${fgbio} DemuxFastqs \
	--inputs ${fastq_dir}/${flow_cell}/${pre_demux_fastq_R1} \
			${fastq_dir}/${flow_cell}/${pre_demux_fastq_R2} \
	--metadata ${sample_sheet} \
	--read-structures 8B12M+T 8M+T \
	--output-type=Fastq \
	--threads ${ncpu} \
	--output ${demux_dir}/fastq \
	--metrics ${demux_dir}/metrics/${metrics_name}_demux_barcode_metrics.txt

# while loop to ensure all steps in this process are done before entering next process
while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
	sleep 60
done
END=$(date +%s)
echo "Demultiplex, time elapsed: $(( $END - $START )) seconds"
echo ""

sed 1d ${sample_sheet} | while read sample_metadata
do

	echo "----------------------------------------------------------------------"
	echo "-------------------  Adapter Trimming with BBduk ---------------------"
	echo "----------------------------------------------------------------------"
	echo ""
	START=$(date +%s)
	Sample_ID=$(echo ${sample_metadata} | tail -n 1 | cut -d ',' -f1)
	pre_trim_fastq_R1=$(ls ${demux_dir}/fastq/${Sample_ID}*_R1.fastq.gz)
	pre_trim_fastq_R2=$(ls ${demux_dir}/fastq/${Sample_ID}*_R2.fastq.gz)

	echo "----------------------------------------------------------------------"
	echo "${BBMap}/bbduk.sh "
	echo "ref=${BBMap}/resources/adapters.fa \ "
	echo "in1=${pre_trim_fastq_R1}"
	echo "in2=${pre_trim_fastq_R2}"
	echo "out1=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz "
	echo "out2=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz "
	echo "ktrim=r k=23 mink=11 hdist=1 trimpolyg=50 tpe tbo "
	echo "----------------------------------------------------------------------"
	echo ""

    #### https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbduk-guide/
	#### /projects/ps-palmer/software/local/src/bbmap/bbduk.sh
	#### ref : removes adapters listed in ${BBMap}/resources/adapters.fa
	#### ktrim=r : right (3') trimming: all adapters and bases to their right are removed
	#### k=23 : contaminant kmer length: removes all reads with a k-length match to seqs in adapters.fa
	#### mink=11 : look for shorter kmers at read tips down to this length
	#### hdist=1 : maximum 'hamming' distance for ref kmers: allows 1 mismatch
	#### trimpolyg=50 : trim poly-G prefixes of at least this length on both ends of reads
	#### tpe : trim pairs evenly: trims both reads to the minimum length of either
	#### tbo : trim by overlap: trims adapters based on where paired reads overlap

	${BBMap}/bbduk.sh \
		ref=${BBMap}/resources/adapters.fa \
		in1=${pre_trim_fastq_R1} \
		in2=${pre_trim_fastq_R2} \
		out1=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz \
		out2=${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz \
		ktrim=r k=23 mink=11 hdist=1 trimpolyg=50 tpe tbo

	while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
		sleep 60
	done
	END=$(date +%s)
	echo "BBduk adapter trimming, time elapsed: $(( $END - $START )) seconds"
	echo ""


	echo "----------------------------------------------------------------------"
	echo "-------------  Quality & Length Trimming with Cutadapt ---------------"
	echo "----------------------------------------------------------------------"
	echo ""
	START=$(date +%s)

	post_trim_fastq_R1=${trimmed_dir}/${Sample_ID}_trimmed_R1.fastq.gz
	post_trim_fastq_R2=${trimmed_dir}/${Sample_ID}_trimmed_R2.fastq.gz

	echo "----------------------------------------------------------------------"
	echo "cutadapt -q 5 "
	echo "--minimum-length 70 --pair-filter=any"
	echo "--cores=0 "
	echo "-o ${post_trim_fastq_R1} "
	echo "-p ${post_trim_fastq_R2}  "
	echo "${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz "
	echo "${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz"
	echo "----------------------------------------------------------------------"
	echo ""

	#### https://cutadapt.readthedocs.io/en/stable/guide.html
	#### q 5 : quality cutoff: trims phred < 5 bases from the 3' end of both (R1,R2) reads
	#### minimum_length 70 : Discard processed reads that are shorter than 70 bp
	#### pair-filter=any : discard read pairs if at least one read (either R1, R2) fulfills the filtering criteria

	source activate hs_rats
	cutadapt -q 5 \
		--minimum-length 70 --pair-filter=any \
		--cores=0 \
		-o ${post_trim_fastq_R1} \
		-p ${post_trim_fastq_R2} \
		${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz
	conda deactivate

	while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
		sleep 60
	done
	END=$(date +%s)
	echo "Cutadapt quality and length trimming, time elapsed: $(( $END - $START )) seconds"
	echo ""

	############################# clean up directory ###############################
	#### clean up directory
	if [ -f ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz ] && [ -f ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz ] && [ -f ${post_trim_fastq_R1} ] && [ -f ${post_trim_fastq_R2} ]; then
		rm ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R1.fastq.gz
		rm ${trimmed_dir}/${Sample_ID}_adapter_trimmed_R2.fastq.gz
		rm -rf ${fastq_dir}
	else
		echo -e "ERROR: something went wrong during demultiplex, or trimming for ${Sample_ID}"
	fi

done

echo ""
echo "----------------------------------------------------------------------"
echo "------------- All demultiplexing and trimming complete! --------------"
echo "----------------------------------------------------------------------"