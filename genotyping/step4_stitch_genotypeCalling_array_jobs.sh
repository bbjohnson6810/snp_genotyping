#!/bin/bash

#### read in declared PBS environment variables
ncpu=${ppn}

#### extract info from argument files
dir_path=$(head -n 9 ${pipeline_arguments} | tail -n 1)
code=$(head -n 11 ${pipeline_arguments} | tail -n 1)
reference_genome=$(head -n 15 ${pipeline_arguments} | tail -n 1)

#### construct more variables based on extracted info
ref_gen=$(echo ${reference_genome} | rev | cut -d '/' -f 1 | cut -d '.' -f 2- | rev)
stitch_path=${dir_path}/${ref_gen}/stitch
code=${code}/genotyping/util

bcftools=$(awk 'BEGIN {count = 0} {if ($1 == "BCFTools") {print $3; exit 0;} else count += 1} END {if (count == NR) {print "ERROR"}}' ${software})
if [ ${bcftools} = "ERROR" ] || [ ! -f ${bcftools} ]; then
    echo "Error: software_location" 
    exit 1
fi

cd $HOME

echo "----------------------------------------------------------------------"
echo "-------------------- HS Rats Genotyping Pipeline ---------------------"
echo "---------------------  Step 4: Genotype Calling   --------------------"
echo "----------------------------------------------------------------------"

echo "----------------------------------------------------------------------"
echo "------------------  Genotype Calling using STITCH   ------------------"
echo "----------------------------------------------------------------------"
START=$(date +%s)

tempdir_chr=${tempdir}

#### each job runs on one chunk in the chunk file
#### number of jobs = number of chunks
#### start the job at the snp position stored on the line in the chunk file
#### that corresponds to the job's array ID 
#### (job 1 starts on line 1, ends on line 2; job 2 starts on line 2, ends on 3; etc)
start_line=${PBS_ARRAYID}
((end_line=start_line+1))
regionStart=$(head -n ${start_line} ${chunk_file} | tail -n 1)
regionEnd=$(head -n ${end_line} ${chunk_file} | tail -n 1)

echo "----------------------------------------------------------------------"
echo "Rscript: ${code}/STITCH.R "
echo "chromosome: ${chr}"
echo "haplotypes (K): ${k}"
echo "nGen: ${nGen}"
echo "nIterations: ${niterations}"
echo "method: ${method}"
echo "stitch path: ${stitch_path}"
echo "temp directory: ${tempdir_chr}"
echo "bamlist: ${bamlist}"
echo "samplenames file: ${sampleNames_file}"
echo "position file: ${posfile}"
echo "regionStart: ${regionStart}"
echo "regionEnd: ${regionEnd}"
echo "nCores: ${nCore}"
echo "reference panels: ${reference_panels}"
echo "----------------------------------------------------------------------"

#### https://github.com/rwdavies/STITCH
#### https://github.com/rwdavies/STITCH/find/master
#### https://github.com/rwdavies/STITCH/blob/master/STITCH/R/functions.R
#### note: not all these arguments are STITCH-specific parameters. The distinction is marked below
#### args[1]: chr = the current chromosome on which the analysis is running
#### args[2]: k = the number of founder haplotypes (8 for HS rats)
#### args[3]: nGen = number of generations since the population was founded (~100 for HS rats)
#### args[4]: niterations = number of STITCH EM iterations
#### args[5]: method = biological model to run
#### args[6]: stitch_path = output directory (stitch: outputdir)
#### args[7]: tempdir_chr = directory for temporary intermediate files (stitch: tempdir) 
#### args[8]: bamlist = path to list of bam file locations
#### args[9]: sampleNames_file = path to list of sample names, in same order as the bamlist
#### args[10]: posfile = path to position file
#### args[11]: regionStart = when running imputation, where to start from
#### args[12]: regionEnd = when running imputation, where to stop
#### args[13]: nCore = number of cores to use
#### args[14]: reference_panels = path to the reference panel of snps on the current chromosome
####           - used to point to two input files:
####           - reference_haplotype_file = path to reference haplotype file in IMPUTE format
####           - (one row per snp, one column per ref haplotype e.g. 1 1 0 0 1 1 1 1 1 1 1 1 1 1 1 1)
####           - reference_legend_file = path to reference haplotype legend file in IMPUTE format
####           - (one row per snp: chr, pos, ref allele, alt allele, e.g. chr1 22585 C A)

source activate hs_rats
Rscript ${code}/STITCH.R \
    ${chr} \
    ${k} \
    ${nGen} \
    ${niterations} \
    ${method} \
    ${stitch_path} \
    ${tempdir_chr} \
    ${bamlist} \
    ${sampleNames_file} \
    ${posfile} \
    ${regionStart} \
    ${regionEnd} \
    ${nCore} \
    ${reference_panels}
conda deactivate

((region_start=regionStart+1))

${bcftools} index -t --threads ${ppn} ${stitch_path}/stitch.${chr}.${region_start}.${regionEnd}.vcf.gz

while [ "$(jobs -rp | wc -l)" -gt 0 ]; do
   sleep 60
done

END=$(date +%s)

echo "Genotype Calling using STITCH, time elapsed: $(( $END - $START )) seconds"
echo ""
echo "----------------------------------------------------------------------"
echo "-------------------  STITCH genotyping complete!   -------------------"
echo "----------------------------------------------------------------------"