#!/bin/bash
#SBATCH -J qc3_final_files
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -t 30:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 30
#SBATCH --mem-per-cpu 16G
#SBATCH -o /tscc/nfs/home/bbjohnson/qc3_final_files-%j.o
#SBATCH -e /tscc/nfs/home/bbjohnson/qc3_final_files-%j.e
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

ncpu=${ppn}

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
util_code=${code}/quality_control/util
metadata=${inputs}/${geno_round}_metadata.csv
genotypes_dir=${dir_path}/genotypes
stitch_path=${dir_path}/${ref_gen}/stitch
stitch_vcf=${stitch_path}/stitch_nonHD_INFO_0.9_rm_discordant.vcf.gz
results_dir=${dir_path}/results
mapping_data=${dir_path}/${ref_gen}/results/mapping_result
mapping_result=${results_dir}/mapping_result
mapping_stats_all=${mapping_data}/${geno_round}_mapping_stats_all.csv
mapping_stats_new=${mapping_data}/mapping_stats_new.csv
genotype_data=${dir_path}/${ref_gen}/results/genotype_result
genotype_result=${results_dir}/genotype_result
previous_flow_cells_metadata=${inputs}/${geno_round}_previous_flowcells_metadata
sample_ids=${inputs}/${geno_round}_sample_ids
chr_names=${inputs}/rn7_2_chr_names
chr_nc=${inputs}/rn7_2_chr_nums


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
echo "-----------    Quality Control 3: Mapping & Genotyping QC    ------------"
echo "-------------------------------------------------------------------------"
echo ""
BEGIN=$(date +%s)

#### get the file prefix for the most recently produced genotype dataset
in_prefix=$(ls ${genotype_data}/*.log  | tail -n 1 | rev |  cut -d '.' -f 2 | cut -d '/' -f 1 | rev)

#### set the output file prefix
n_bams=$(cat ${bamlist} | wc -l)
datestamp=$(date +%y%m%d)
out_prefix=hs_rats_${geno_round}_n${n_bams}_20${datestamp}

#### create a sample_id list in the same order they occur in the .vcf files
#### (to be used while making the genotype log to produce sample subsets)
## ${bcftools} query -l ${stitch_path}/stitch_HD.vcf.gz > ${results_dir}/stitch_HD_sample_ids
## ${bcftools} query -l ${genotypes_dir}/round10_1_raw_nofilters.vcf.gz > ${results_dir}/stitch_HD_sample_ids

echo "------------------ Genotyping log and quality control -------------------"
echo ""
START=$(date +%s)
source activate hs_rats


#### Execute quality control tests to produce a genotyping log
python3 ${util_code}/make_genotype_log.py \
	-i ${metadata} \
	-m ${mapping_stats_all} \
	-x ${genotype_data}/${in_prefix}.smiss \
	-z ${genotype_data}/${in_prefix}.het \
	-s ${sample_ids} \
	-r ${geno_round} \
	-o ${results_dir}/${out_prefix}

gtype_log_produced=$?

save the genotype log
gtype_log=${results_dir}/${out_prefix}_genotype_log.csv
echo "gtype_log: ${gtype_log}"
# conda deactivate
END=$(date +%s)

if [ ${gtype_log_produced} == 0 ]; then
	echo ""
	echo "QC results saved to ${results_dir}/${out_prefix}_genotype_log.csv"
	echo ""
else 
	echo ""
	echo "Error conducting QC (see traceback above). No genotype log produced!"
	echo ""
fi
echo "Genotype log time elapsed: $(( $END - $START )) seconds"
echo ""


echo "---------------------- Save QC statistics ---------------------"
echo ""
START=$(date +%s)
source activate hs_rats
ls ${util_code}/summary_stats.py
python3 ${util_code}/summary_stats.py \
	-g ${gtype_log} \
	-i ${genotype_data}/${geno_round}_INFO_all_snps \
	-a ${genotype_data}/${in_prefix}.afreq \
	-w ${genotype_data}/${in_prefix}.hardy \
	-x ${genotype_data}/${in_prefix}.vmiss \
	-r ${geno_round} \
	-o ${results_dir}/${geno_round}

conda deactivate
END=$(date +%s)
echo "QC summary stats time elapsed: $(( $END - $START )) seconds"
echo ""


echo "---------------------- Plot QC statistics ---------------------"
echo ""
START=$(date +%s)
source activate hs_rats

echo "Plotting read mapping..."
python3 ${util_code}/plot_reads_mapped.py \
	-i ${mapping_stats_new} \
	-o ${mapping_result}/${geno_round}_mapstats

echo "Plotting reads mapped by chromosome..."
python3 ${util_code}/plot_reads_mapped_chr.py \
	-i ${mapping_stats_new} \
	-s ${metadata} \
	-o ${mapping_result}/${geno_round}_mapstats_chr

echo "Plotting general mapping statistics..."
python3 ${util_code}/plot_general_mapstats.py \
	-i ${gtype_log} \
	-o ${mapping_result}/${geno_round}_mapstats

echo "Plotting sex QC: read mapping on chrX vs chrY..."
python3 ${util_code}/plot_sexQC_chrX_vs_chrY.py \
	-i ${gtype_log} \
	-o ${mapping_result}/${geno_round}_qc_sex

echo "Plotting INFO score histograms..."
python3 ${util_code}/plot_INFO_hist.py \
	-i ${genotype_data}/${geno_round}_INFO_all_snps \
	-o ${genotype_result}/${geno_round}_INFO_all

python3 ${util_code}/plot_INFO_hist.py \
	-i ${genotype_data}/${geno_round}_INFO_0.9_snps \
	-o ${genotype_result}/${geno_round}_INFO_0.9

echo "Plotting INFO score chromosome maps"
Rscript ${util_code}/plot_info_chrs.R \
	${genotype_data}/${geno_round}_INFO_all_snps \
	${genotype_data}/${geno_round}_INFO_0.9_snps \
	${genotype_result}/${geno_round}

### to do: plot snp density map like info scores ###
# look at genotypes_SNPs_density.py

###########################################################
## convert this code to accommodate new mapstats file ###
echo "Plotting "
python3 ${util_code}/genotypes_missing_vs_mapped_reads.py \
	-r ${mapping_stats_all} \
	-m ${genotype_data}/${in_prefix}.smiss \
	-o ${genotype_result}/${geno_round}
############################################################

echo "Plotting sample heterozygosity vs. missingness"
## be sure to set upper and lower heterozygosity cutoffs for plots as needed
python3 ${util_code}/plot_sample_het_vs_missing.py \
	-i ${gtype_log} \
	-m "ddGBS" \
	-u 4 \
	-l 4 \
	-o ${genotype_result}/${geno_round}_sample_qc_het_vs_missing_ddGBS

python3 ${util_code}/plot_sample_het_vs_missing.py \
	-i ${gtype_log} \
	-m "lcWGS" \
	-u 5 \
	-l 4 \
	-o ${genotype_result}/${geno_round}_sample_qc_het_vs_missing_lcWGS

echo "Plotting SNP statistics"
python3 ${util_code}/plot_snp_stats.py \
	-i ${genotype_data}/${geno_round}_INFO_0.9_snps \
	-a ${genotype_data}/${in_prefix}.afreq \
	-w ${genotype_data}/${in_prefix}.hardy \
	-x ${genotype_data}/${in_prefix}.vmiss \
	-o ${genotype_result}/${geno_round}_snpstats_INFO_0.9

python3 ${util_code}/plot_snp_stats.py \
	-i ${genotype_data}/${geno_round}_INFO_all_snps \
	-a ${genotype_data}/${in_prefix}.afreq \
	-w ${genotype_data}/${in_prefix}.hardy \
	-x ${genotype_data}/${in_prefix}.vmiss \
	-o ${genotype_result}/${geno_round}_snpstats_INFO_all

echo "Plotting genotypes PCA"
python3 ${util_code}/pca_raw_genotypes.py \
	-c ${genotype_data}/${in_prefix}.eigenvec \
	-v ${genotype_data}/${in_prefix}.eigenval \
	-m ${metadata} \
	-o ${genotype_result}/pca/raw/${geno_round}_raw

conda deactivate

END=$(date +%s)
echo "QC plotting time elapsed: $(( $END - $START )) seconds"
echo ""




echo ""
echo "-------------------------------------------------------------------------"
echo "------------------------ Genotyping QC complete! ------------------------"
echo "-------------------------------------------------------------------------"
echo ""

