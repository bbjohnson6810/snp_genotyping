#!/usr/win/env python3

# This script conducts quality control checks on read mapping and sample genotypes,
# and produces the final genotype report and quality control summary

# Usage: python3 make_genotype_log.py -i sample_sheet -m mapping_stats -x missingness_file \
#                                     -z het_file -s stitch_sample_ids -o output_file_prefix

import os
import sys
import numpy as np
import pandas as pd
import datetime
from optparse import OptionParser

print('working')

def help():
	print("====== make_genotype_log.py =====")
	print("Produce a genotyping log with summary stats and QC results")
	print("-i <sample sheet>                    the sample sheet with metadata for all genotyped samples")
	print("-m <mapping statistics file>         the csv file with read mapping statistics for all genotyped samples")
	print("-x <plink sample missingness file>   the plink .smiss file with sample-based data on missing genotypes")
	print("-z <plink heterozygosity file>       the plink .het file with sample-based heterozygosity data")
	print("-s <sample ID file>                  file listing sample_ids in the stitch_HD.vcf file, in the same order as the vcf file")
	print("-o <output file prefix>              the output file prefix: same as the <plink_prefix> used for genotype data")
	print("Usage: python3 make_genotype_log.py -i sample_sheet_file -m mapping_stats_file -o output_file_prefix")
	sys.exit()

#### EDIT ####
if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-H', action="store_true", dest="show_help", help="show this help message and exit")
    parser.add_option('-i', type="string", nargs=1, dest="sample_sheet", help="<sample sheet for all genotyped rats>")
    parser.add_option('-m', type="string", nargs=1, dest="map_stats", help="<mapping statistics file>")
    parser.add_option('-x', type="string", nargs=1, dest="missingness_file", help="<plink .smiss file>")
    parser.add_option('-z', type="string", nargs=1, dest="het_file", help="<plink .het file>")
    parser.add_option('-s', type="string", nargs=1, dest="stitch_sample_ids", help="<STITCH sample IDs>")
    parser.add_option('-r', type="string", nargs=1, dest="gtyping_round", help="<genotyping pipeline round number>")
    parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
    options, args = parser.parse_args()
        
        
    if len(sys.argv) == 1 or options.show_help:
        help()
    if options.out_file_prefix != None:
        output_file_prefix = options.out_file_prefix
    else:
        raise "Please provide an output file prefix"
    if options.sample_sheet != None:
        sample_sheet = options.sample_sheet
    else:
        raise "Please provide a sample sheet file"
    if options.map_stats != None:
        map_stats = options.map_stats
    else:
        raise "Please provide a mapping statistics file"
    if options.missingness_file != None:
        missingness_file = options.missingness_file
    else:
        raise "Please provide a .smiss file with genotype missingness data"
    if options.het_file != None:
        het_file = options.het_file
    else:
        raise "Please provide a .het file with heterozygosity data"
    if options.stitch_sample_ids != None:
            stitch_ids = options.stitch_sample_ids
    else:
        raise "Please provide a list of STITCH sample IDs"
    if options.gtyping_round != None:
        gtype_round = options.gtyping_round
    else:
        raise "Please provide the current genotyping round"


    # read in metadata
    md = pd.read_csv(sample_sheet)

    # create a new column for sequencing method
    md['seq_method'] = np.where(md['library_name'].str.contains('GBS'), 'ddGBS', 'lcWGS')

    # rearrange columns 
    md = md[['sample_id', 'library_name', 'rfid', 'barcode', 'pcr_barcode', 'project_name', 'runid', 'seq_method','organism', 'strain', 
             'coatcolor', 'sex', 'dams', 'sires']]
    md.rename(columns={'Sample_ID':'sample_id', 'library_name':'library', 'dams':'dam', 'sires':'sire'}, inplace=True)

    # read mapping statistics
    mapstats = pd.read_csv(map_stats)

    # save columns to merge to the metadata file
    map_cols = ['sample_id', 'total_reads', 'primary_reads', 'total_mapped', 'primary_mapped', 'duplicate_primary_reads', 
                    'chr1_mapped', 'chr2_mapped', 'chr3_mapped', 'chr4_mapped', 'chr5_mapped', 'chr6_mapped', 'chr7_mapped', 
                    'chr8_mapped', 'chr9_mapped', 'chr10_mapped', 'chr11_mapped', 'chr12_mapped', 'chr13_mapped', 
                    'chr14_mapped', 'chr15_mapped', 'chr16_mapped', 'chr17_mapped', 'chr18_mapped','chr19_mapped', 
                    'chr20_mapped', 'chrX_mapped', 'chrY_mapped', 'chrM_mapped']
    mapstats = mapstats[map_cols]

    # merge mapping statistics to the metadata file
    md = pd.merge(md, mapstats, how='left', on='sample_id').reset_index(drop=True)

    # drop sample IDs from the list of mapstats columns
    map_cols = map_cols[1:]

    # transform all read counts for easier interpretability: divide by 1M
    for col in md[map_cols]:
        md[col] = md[col] / 1e6
    

    #### READ MAPPING QC ####
    
    # set a threshold of 2M reads per sample
    # samples with < 2M reads fail read mapping QC
    mapped_reads_threshold = 2
    md['qc_mapped_reads'] = md['primary_mapped'].apply(lambda x: 
                                                       'pass' if x >= mapped_reads_threshold else 'fail')

    # calculate the percentage of reads that mapped
    md['pct_mapped'] = md['primary_mapped'] / md['primary_reads'] * 100

    # calculate the percentage of duplicate reads
    md['pct_duplicated'] = md['duplicate_primary_reads'] / md['primary_reads'] * 100


    #### SEX QC ####

    # QC function: sets cutoffs for read percentages    
    def sex_QC(x, y, sex):
        if isinstance(sex, str):  # check if 'sex' is a string
            sex = sex.strip()  # remove leading and trailing whitespaces
            if not sex:  # check if 'sex' is empty after stripping
                return 'not_tested'
            elif sex == 'F':
                if y >= x/20.5 + 0.07 + 0.02:
                    return 'fail'
                elif y >= x/20.5 + 0.07 - 0.02:
                    return 'suspect'
                else:
                    return 'pass'
            elif sex == 'M':
                if y <= x/20.5 + 0.07 - 0.02:
                    return 'fail'
                elif y <= x/20.5 + 0.07 + 0.02:
                    return 'suspect'
                else:
                    return 'pass'
        else:
            return 'not_tested'

    # calculate percentage of reads mapped to either sex chromosome
    md['pct_chrX'] = md['chrX_mapped'] / md['total_mapped'] * 100
    md['pct_chrY'] = md['chrY_mapped'] / md['total_mapped'] * 100

    # sex QC designations
    md['qc_sex'] = md.apply(lambda row: sex_QC(row['pct_chrX'], row['pct_chrY'], row['sex']), axis=1)


    #### MISSING RATE QC ####
    
    # add missingness data to metadata
    sample_missing = pd.read_csv(missingness_file, sep='\t', usecols=['IID', 'F_MISS'])
    sample_missing = sample_missing.rename(columns={'IID': 'sample_id', 'F_MISS': 'missing_rate'})
    md = pd.merge(md, sample_missing, how='left', on='sample_id').reset_index(drop=True)

    # missingness QC: fail if a sample has >10% of snps missing
    missing_rate_threshold = 0.1
    md['qc_missing_rate'] = md['missing_rate'].apply(lambda x:
                                                        'pass' if x <= missing_rate_threshold else 'fail')


    #### HETEROZYGOSITY QC ####

    sample_het = pd.read_csv(het_file, sep='\t', usecols=['IID', 'O(HOM)', 'OBS_CT'])
    sample_het['O(HOM)'] = pd.to_numeric(sample_het['O(HOM)'])
    sample_het['OBS_CT'] = pd.to_numeric(sample_het['OBS_CT'])
 
    # heterozygosity = (# of snp genotypes - # of homozygous genotypes) / # of snp genotypes
    #                = # heterozygous snps / # of snps = frequency of snps that are heterozygous
    sample_het['heterozygosity'] = (sample_het['OBS_CT'] - sample_het['O(HOM)'])/sample_het['OBS_CT']
    sample_het = sample_het.rename(columns={'IID': 'sample_id'})
        
    md = pd.merge(md, sample_het[['sample_id', 'heterozygosity']], how='left', on='sample_id').reset_index(drop=True)

    # heterozygosity QC:
    # extreme heterozygosity rates (> 4std from the mean) fail if they have low read mapping
    # samples with extreme heterozygosity but passing read mapping QC are flagged
    # samples with 'normal' heterozygosity but low read mapping are flagged

    founder_ids = ['ACI', 'BN', 'BUF', 'F344', 'M520', 'MR', 'WKY', 'WN', 
                'M520N', 'WKYN']

    # heterozygosity summary stats based on sequencing method
    het_mean_GBS = md[md['seq_method']=='ddGBS']['heterozygosity'].mean()
    het_std_GBS = md[md['seq_method']=='ddGBS']['heterozygosity'].std()
    # het_mean_lcWGS = md[md['seq_method']=='lcWGS']['heterozygosity'].mean()
    # het_std_lcWGS = md[md['seq_method']=='lcWGS']['heterozygosity'].std()
    het_mean_lcWGS = md[(md['seq_method']=='lcWGS') & (~md['rfid'].isin(founder_ids))]['heterozygosity'].mean()
    het_std_lcWGS = md[(md['seq_method']=='lcWGS') & (~md['rfid'].isin(founder_ids))]['heterozygosity'].std()

    def qc_heterozygosity(row):
        
        if row['seq_method']=='ddGBS':
            het_mean = het_mean_GBS
            het_std = het_std_GBS
        else:
            het_mean = het_mean_lcWGS
            het_std = het_std_lcWGS
        
        # heterozygosity rates > 4std from the mean fail if they have low read mapping
        if row['heterozygosity'] > het_mean+4*het_std and row['missing_rate'] > missing_rate_threshold:
            return 'fail'
        elif row['heterozygosity'] < het_mean-4*het_std and row['missing_rate'] > missing_rate_threshold:
            return 'fail'

        # founder genotypes with low heterozygosity can pass
        elif (row['heterozygosity'] < het_mean - 4 * het_std) and (row['rfid'] in founder_ids):
            return 'pass'

        # heterozygosity rates > 4std from the mean are flagged if they have good read mapping
        elif row['heterozygosity'] > het_mean+4*het_std and row['missing_rate'] <= missing_rate_threshold:
            return 'suspect'
        elif row['heterozygosity'] < het_mean-4*het_std and row['missing_rate'] <= missing_rate_threshold:
            return 'suspect'

        # all other samples pass QC
        else:
            return 'pass'
    

    # run heterozygosity QC
    md['qc_heterozygosity'] = md.apply(lambda row: qc_heterozygosity(row), axis=1)


    #### GENOTYPE LOG ####
    # track genotyping counts and QC performance for each sample

    # flag when a sample has passed all QC steps
    def all_qc(row):
        
        if (row['qc_sex'] == 'pass' or row['qc_sex'] == 'not_tested') \
            and row['qc_missing_rate'] == 'pass' and row['qc_heterozygosity'] == 'pass':
            return 'pass'
        else: 
            return 'fail'

    md['qc_all'] = md.apply(all_qc, axis=1)

    # flag if a given sample has ever been re-genotyped
    all_rfids = md['rfid'].tolist()
    dup_samples = set([x for x in all_rfids if all_rfids.count(x) > 1])
    md['duplicate_rfid'] = md['rfid'].isin(dup_samples).astype(int)

    # count the number of times each rfid has been genotyped
    redo_counts = md['rfid'].value_counts()
    md['times_genotyped'] = md['rfid'].apply(lambda row: redo_counts[row]).astype(int)

    # count the number of times each rfid has passed QC
    pass_counts = md[md['qc_all'] == 'pass'].groupby('rfid', as_index=False).size()
    pass_counts = pass_counts.rename(columns={"size": "times_passed"})
    md = md.merge(pass_counts, how = 'left', on='rfid')
    md['times_passed'] = md['times_passed'].fillna(0).astype(int)

    # count the number of times each rfid has failed QC
    fail_counts = md[md['qc_all'] == 'fail'].groupby('rfid', as_index=False).size()
    fail_counts = fail_counts.rename(columns={"size": "times_failed"})
    md = md.merge(fail_counts, how = 'left', on='rfid')
    md['times_failed'] = md['times_failed'].fillna(0).astype(int)

    # keep only samples that passed sex, missingness, and heterozygosity QC
    qc_all_pass = md[md['qc_all']=='pass'].reset_index(drop=True)

    # remove duplicates in the samples that passed QC 
    # then reorder metadata by sample id
    qc_all_pass_samples = qc_all_pass['rfid'].tolist()
    dup_samples = set(qc_all_pass[qc_all_pass['times_passed']>1]['rfid'])
    dup_md = qc_all_pass[qc_all_pass['rfid'].isin(
        dup_samples)].sort_values(by=['rfid', 'library']).reset_index(drop=True)

    # keep sample runs with the most mapped reads
    dup_keep =[]
    for s in dup_samples:
        temp_sample_id = dup_md[dup_md['rfid'] == s]['sample_id'].tolist()
        temp_mapped_reads = dup_md[dup_md['rfid'] == s]['primary_mapped'].tolist()
        dup_keep.append(temp_sample_id[temp_mapped_reads.index(max(temp_mapped_reads))])
    dup_drop = dup_md[~dup_md['sample_id'].isin(dup_keep)]['sample_id'].tolist()

    qc_all_pass_no_dup = qc_all_pass.drop(qc_all_pass[qc_all_pass['sample_id'].isin(
        dup_drop)].index).reset_index(drop=True)


    #### SAMPLE USAGE ####

    # function to designate how each sample should be used
    def sample_usage(row):
        
        # save all instances in which this rfid was genotyped
        sample_md = md[md['rfid']==row['rfid']]
        
        # save counts of the number of time the sample passed or failed all QC
        qc_counts = sample_md['qc_all'].value_counts().rename_axis('outcome').reset_index(name='count')
        
        # if an rfid has only ever passed or failed, add a row denoting zero fails or passes
        if set(qc_counts['outcome']) == {'pass'}:
            qc_counts.loc[1] = ['fail',0]
        elif set(qc_counts['outcome']) == {'fail'}:
            qc_counts.loc[1] = ['pass',0]

        times_passed = int(qc_counts[qc_counts['outcome']=='pass']['count'].iloc[0])
        times_failed = int(qc_counts[qc_counts['outcome']=='fail']['count'].iloc[0])

        # save mapped read counts for all instances this rfid was sequenced
        sample_reads = row['primary_mapped']
        if row['times_passed'] > 1:
            reads_set = sample_md[sample_md['qc_all']=='pass']['primary_mapped'].tolist()
        else:
            reads_set = sample_md['primary_mapped'].tolist()
        
        ## sample usage designations ##

        # samples that passed all QC on their 1st and only run are saved for analyses
        if (row['times_genotyped'] == 1) & (row['times_failed'] == 0):
            return 'analysis'
        
        # samples that passed all QC only once (of multiple runs) are saved for analyses
        elif (row['times_genotyped'] > 1) & (times_passed == 1) & (row['qc_all'] == 'pass'):
            return 'analysis'
        
        # for samples that passed all QC multiple times, 
        # only the run with the most mapped reads is saved for analysis
        elif (row['times_genotyped'] > 1) & (times_passed > 1) & (sample_reads == max(reads_set)) \
        & (row['qc_all'] == 'pass'):
            return 'analysis'

        # samples that passed QC multiple times but aren't used for analysis 
        # can be kept in future genotyping rounds (to help with imputation)
        elif (row['times_genotyped'] > 1) & (times_passed > 1) & (sample_reads < max(reads_set)) \
        & (row['qc_all'] == 'pass'):
            return 'keep'
        
        # samples that failed QC on their 1st (and only) run need to be re-sequenced
        elif (row['times_genotyped'] == 1) & (times_passed == 0):
            return 'resequence'
        
        # samples that failed QC 2 out of 2 times need to be re-extracted 
        # using tail tissue and then resequenced using the new extraction
        elif (row['times_genotyped'] == 2) & (times_passed == 0):
            return 'reseq_tail'

        # samples that failed QC 3 out of 3 times need to be dropped from all future genotyping rounds
        elif (row['times_genotyped'] == 3) & (times_passed == 0):
            return 'drop'

        # for samples that were run multiple times and passed QC at least once,
        # all failed runs can be dropped from future genotyping rounds
        elif (row['times_genotyped'] > 1) & (times_passed > 0) & (row['qc_all'] == 'fail'):
            return 'drop'

        else:
            return 'LEFTOVER-CHECK'

    # designate sample usages
    md['sample_use'] = md.apply(lambda x: sample_usage(x), axis=1)


    # add  pipeline version info to metadata
    md['pipeline_version'] = 'Hybrid_v1.0.2'
    md['pipeline_round'] = gtype_round

    # rearrange columns so QC columns are together
    # print(md.columns.tolist()[42])
    out_cols = ['sample_id', 'library', 'rfid', 'project_name', 'barcode', 'pcr_barcode','runid', 'seq_method', 
                'pipeline_version', 'pipeline_round', 'organism', 'strain', 'coatcolor', 'sex', 'dam', 'sire',
                'primary_reads', 'duplicate_primary_reads','primary_mapped', 'pct_mapped', 'pct_chrX', 'pct_chrY', 
                'missing_rate', 'heterozygosity','qc_sex', 'qc_missing_rate', 'qc_heterozygosity', 'qc_all', 
                'duplicate_rfid', 'times_genotyped','times_passed', 'times_failed', 'sample_use']

    md = md[out_cols]


    #### OUTPUT FILES ####

    #### genotype log table for all samples
    md.to_csv(output_file_prefix + '_genotype_log.csv',index=False)

    #### log of all samples that pass final QC for use in analyses
    md_for_analysis = md[md['sample_use']=='analysis']
    md_for_analysis.to_csv(f'{output_file_prefix}_use_for_analysis_n{len(md_for_analysis)}_genotype_log.csv', index=False)

    #### log of all samples that fail final QC and need to be re-sequenced or dropped
    md_to_drop = md[~md['sample_use'].isin(['analysis', 'keep'])].copy()
    md_to_drop.reset_index(drop=True, inplace=True)
    md_to_drop.to_csv(f'{output_file_prefix}_reseq_drop_n{len(md_to_drop)}_genotype_log.csv', index=False)

    # create sample files to subset and rename samples in .vcf files
    # note: all MUST be in the same order as in STITCH .vcf files
    stitch_ids = pd.read_csv(stitch_ids,header=None)
    stitch_ids = stitch_ids[0].tolist()
    md_ids = md['sample_id'].tolist()
    md_use_ids = md_for_analysis['sample_id'].tolist()
    stitch_use_ids = [x for x in stitch_ids if x in md_use_ids]

    # function to reorder sample_ids (list1) based on the order of stitch sample_ids (list2)
    def reorder_list(list1, list2):
        # create a mapping between the elements of list2 and their indices
        mapping = {value: index for index, value in enumerate(list2)}
        
        # sort list1 based on the indices from the mapping
        sorted_list1 = sorted(list1, key=lambda x: mapping.get(x, float('inf')))
        
        return sorted_list1

    # reorder sample_ids
    md_use_ids = reorder_list(md_use_ids, stitch_use_ids)
 
    # reorder rfids corresponding to stitch sample_ids
    rfid_dict = dict(zip(md_for_analysis['sample_id'],md_for_analysis['rfid']))
    md_use_rfids = [rfid_dict[key] for key in stitch_use_ids]

    #### list of sample IDs - used by Bcftools view -S to subset samples
    # (to remove QC-failed samples from the VCF file)
    pd.DataFrame(md_use_ids).to_csv(f'{output_file_prefix}_use_for_analysis_n{len(md_for_analysis)}_sample_ids', 
                    index=False, sep=' ', header=False)

    #### lists of sample IDs and rfids - used by Bcftools reheader -S to rename samples
    # (to remove the library ID in the sample names of the VCF file)
    pd.DataFrame({'sample_id':md_use_ids, 'rfid':md_use_rfids}).to_csv(
        f'{output_file_prefix}_use_for_analysis_n{len(md_for_analysis)}_sample_rename',index=False, sep=' ', header=False)
    pd.DataFrame({'rfid':md_use_rfids}).to_csv(
        f'{output_file_prefix}_use_for_analysis_n{len(md_for_analysis)}_rfids',index=False, sep=' ', header=False)
