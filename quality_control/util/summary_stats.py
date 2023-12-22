#!/usr/win/env python3

# This script prints genotyping summary statistics to a text file

# Usage: python3 make_genotype_log.py -i sample_sheet -m mapping_stats -x missingness_file \
#                                     -z het_file -s stitch_sample_ids -o output_file_prefix

import os
import sys
import numpy as np
import pandas as pd
import datetime
from optparse import OptionParser


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

# function to produce a dataset of snp stats
def make_snp_df(info_file, afreq_file, hwe_file, vmiss_file):
    
    print('Setting up SNP data...')

    # read in various snp data
    info = pd.read_csv(info_file, delimiter='\t', dtype=str, usecols=['CHROM', 'POS', 'INFO_SCORE'])
    afreq = pd.read_csv(afreq_file, delimiter='\t', dtype=str, usecols=['ID', 'ALT_FREQS'])
    hwe = pd.read_csv(hwe_file, delimiter='\t', dtype=str, usecols=['ID', 'P'])
    vmiss = pd.read_csv(vmiss_file, delimiter='\t', dtype=str, usecols=['ID', 'F_MISS'])

    # clean up snp data
    info['CHROM'] = info['CHROM'].apply(lambda x: x.split('chr')[1])
    info['INFO_SCORE'] = pd.to_numeric(info['INFO_SCORE'])
    info['POS'] = pd.to_numeric(info['POS'])
    # info = info.map(str)
    print('info:')
    print(info.head())

    afreq['CHROM'] = afreq['ID'].apply(lambda x: x.split(':')[0])
    afreq['POS'] = afreq['ID'].apply(lambda x: x.split(':')[1])
    afreq['POS'] = pd.to_numeric(afreq['POS'])
    afreq['ALT_FREQS'] = pd.to_numeric(afreq['ALT_FREQS'])
    afreq['MAF'] = afreq['ALT_FREQS'].apply(lambda x: x if x <=0.5 else 1-x)
    # afreq = afreq.map(str)
    print('afreq:')
    print(afreq.head())

    hwe['CHROM'] = hwe['ID'].apply(lambda x: x.split(':')[0])
    hwe['POS'] = hwe['ID'].apply(lambda x: x.split(':')[1])
    hwe['POS'] = pd.to_numeric(hwe['POS'])
    hwe['P'] = pd.to_numeric(hwe['P'])
    hwe['-log(P)'] = -np.log10(hwe['P'])
    print('hwe:')
    print(hwe.head())

    vmiss['CHROM'] = vmiss['ID'].apply(lambda x: x.split(':')[0])
    vmiss['POS'] = vmiss['ID'].apply(lambda x: x.split(':')[1])
    vmiss['POS'] = pd.to_numeric(vmiss['POS'])
    vmiss['F_MISS'] = pd.to_numeric(vmiss['F_MISS'])
    vmiss = vmiss.rename(columns={'F_MISS':'missing_rate'})
    print('vmiss:')
    print(vmiss.head())

    # # merge all data
    snp_dat = pd.merge(info, afreq[['CHROM', 'POS', 'MAF']], how='left', 
                       on=['CHROM', 'POS']).reset_index(drop=True)
    snp_dat = pd.merge(snp_dat, hwe[['CHROM', 'POS', '-log(P)']], 
                       how='left', on=['CHROM', 'POS']).reset_index(drop=True)
    snp_dat = pd.merge(snp_dat, vmiss[['CHROM', 'POS', 'missing_rate']], 
                       how='left', on=['CHROM', 'POS']).reset_index(drop=True)
    snp_dat = snp_dat.sort_values(by=['CHROM', 'POS']).reset_index(drop=True)

    print('SNP dataset complete')
    print(snp_dat.head())
    return snp_dat

#### EDIT ####
if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-H', action="store_true", dest="show_help", help="show this help message and exit")
    parser.add_option('-g', type="string", nargs=1, dest="gtype_log", help="<current genotyping log>")
    parser.add_option('-r', type="string", nargs=1, dest="gtype_round", help="<genotyping pipeline round number>")
    parser.add_option('-i', type="string", nargs=1, dest="info_scores", help="<INFO scores for all SNPs>")
    parser.add_option('-a', type="string", nargs=1, dest="afreq_file", help="<plink .afreq file>")
    parser.add_option('-w', type="string", nargs=1, dest="hwe_file", help="<plink .hardy file>")
    parser.add_option('-x', type="string", nargs=1, dest="vmiss_file", help="<plink .vmiss file")
    parser.add_option('-o', type="string", nargs=1, dest="output_prefix", help="<output file prefix>")
    options, args = parser.parse_args()
        
        
    if len(sys.argv) == 1 or options.show_help:
        help()
    if options.gtype_log != None:
        gtype_log = options.gtype_log
    else:
        raise "Please provide the current genotyping log"
    if options.gtype_round != None:
        gtype_round = options.gtype_round
    else:
        raise "Please provide the current genotyping round"
    if options.info_scores != None:
        info_all = options.info_scores
    else:
        raise "Please provide INFO scores for all SNPs"
    if options.output_prefix != None:
        output_prefix = options.output_prefix
    else:
        raise "Please provide an output file prefix"
    if options.afreq_file != None:
        afreq = options.afreq_file
    else:
        raise "Please provide a plink .afreq file"
    if options.hwe_file != None:
        hwe = options.hwe_file
    else:
        raise "Please provide a plink .hardy file"
    if options.vmiss_file != None:
        vmiss = options.vmiss_file
    else:
        raise "Please provide a plink .vmiss file"


    # read in metadata
    md = pd.read_csv(gtype_log)

    # create snp dataframe
    snp_df = make_snp_df(info_all, afreq, hwe, vmiss)
    snps_use = snp_df[snp_df['INFO_SCORE'] > 0.9]
    snps_maf = snps_use[snps_use['MAF'] > 0.01]
    snps_miss = snps_use[snps_use['missing_rate'] < 0.1]
    snps_hwe = snps_use[snps_use['-log(P)'] < 6]
    snps_filtered = pd.merge(snps_maf, snps_miss, how='left', 
                       on=['CHROM', 'POS']).reset_index(drop=True)
    snps_filtered = pd.merge(snps_filtered, snps_hwe, how='left', 
                       on=['CHROM', 'POS']).reset_index(drop=True)

    info_all = pd.read_csv(info_all)

    #### QC SUMMARY ####

    # file to output summary stats
    summary_file = f'{output_prefix}_summary_stats.txt'
    ct = datetime.datetime.now()

    # write various summary stats to file
    with open(summary_file, 'w') as outfile:
        
        print('#### GENOTYPING SUMMARY STATS ####', file=outfile)
        print('', file=outfile)
        print(f'Summary produced on {ct}', file=outfile)
        print('', file=outfile)
        
        # sample counts
        print('', file=outfile)
        print('## SAMPLE COUNTS ##', file=outfile)
        print('', file=outfile)
        print(f'Total number of samples: {len(md)}', file=outfile)
        print(f'    Total male samples: {len(md[md["sex"] == "M"])}', file=outfile)
        print(f'    Total female samples: {len(md[md["sex"] == "F"])}', file=outfile)
        print('', file=outfile)
        print(f'Total ddGBS samples: {len(md[md["seq_method"]=="ddGBS"])}', file=outfile)
        print(f'    Male ddGBS: {len(md[(md["seq_method"]=="ddGBS") & (md["sex"]=="M")])}', file=outfile)
        print(f'    Female ddGBS: {len(md[(md["seq_method"]=="ddGBS") & (md["sex"]=="F")])}', file=outfile)
        print('', file=outfile)
        print(f'Total lcWGS samples: {len(md[md["seq_method"]=="lcWGS"])}', file=outfile)
        print(f'    Male lcWGS: {len(md[(md["seq_method"]=="lcWGS") & (md["sex"]=="M")])}', file=outfile)
        print(f'    Female lcWGS: {len(md[(md["seq_method"]=="lcWGS") & (md["sex"]=="F")])}', file=outfile)
        print('', file=outfile)
        print(f'Total number of unique samples: {len(md["rfid"].unique())}', file=outfile)
        print('', file=outfile)

        # final QC
        print('', file=outfile)
        print('## OVERALL QC ##', file=outfile)
        print('', file=outfile)
        print(f'Total number of samples passing QC: {len(md[md["qc_all"]=="pass"])} ({len(md[md["qc_all"]=="pass"])/len(md)*100:.2f} %)', file=outfile)
        print('', file=outfile)
        all_samples = md['rfid'].tolist()
        redone_times = set([all_samples.count(x) for x in all_samples])
        redone_samples={}
        for i in redone_times:
            redone_samples[i] = set([x for x in all_samples if all_samples.count(x) == i])
            print('Number of samples that have been redone ' + str(i-1) + ' time(s): '
            + str(len(redone_samples[i])), file=outfile)
        print('', file=outfile)
        print('Total number of samples that passed sex, missingness, and heterozygosity QCs (including redos): ' + str(len(md[md["qc_all"]=="pass"])), file=outfile)
        nodups = list(set(md[md["qc_all"]=="pass"]))
        print('Number of unique samples that passed sex, missingness, and heterozygosity QCs: ' + str(len(nodups)), file=outfile)
        pass_all_qc_samples = md[md['qc_all']=='pass']['rfid'].tolist()
        failed_samples = {}
        for i in redone_times:
            failed_samples[i] = [x for x in redone_samples[i] if x not in pass_all_qc_samples]
            print('Number of samples that have been redone ' + str(i-1) + ' time(s) but still failed sex, missing or het QC: '
                +str(len(failed_samples[i])), file=outfile)
        print('', file=outfile)
        
        # read mapping
        mapped_reads_threshold = 2
        print('', file=outfile)
        print('## READ MAPPING ##', file=outfile)
        print('', file=outfile)
        print(f'QC mapped reads threshold: {mapped_reads_threshold} million', file=outfile)
        print(f'    # of passing samples: {len(md[md["qc_mapped_reads"] == "pass"])}', file=outfile)
        print(f'    # of failed samples: {len(md[md["qc_mapped_reads"] == "fail"])}', file=outfile)
        print('', file=outfile)

        # read counts
        print('# read counts #', file=outfile)
        print('', file=outfile)
        
        for sm in md['seq_method'].unique():
            print(f'Number of mapped reads (million), {sm} samples:', file=outfile)
            print(f'    Mean: {md[md["seq_method"]==sm]["primary_mapped"].mean():.2f}', file=outfile)
            print(f'    SD: {md[md["seq_method"]==sm]["primary_mapped"].std():2f}', file=outfile)
            print('', file=outfile)

        # mapping rates
        print('# mapping rates #', file= outfile)
        print('', file=outfile)
        
        for sm in md['seq_method'].unique():
            print(f'Mapping rate (%), {sm} samples:', file=outfile)
            print(f'    Mean: {md[md["seq_method"]==sm]["pct_mapped"].mean():.2f}', file=outfile)
            print(f'    Median: {md[md["seq_method"]==sm]["pct_mapped"].median():.2f}', file=outfile)
            print(f'    SD: {md[md["seq_method"]==sm]["pct_mapped"].std():.2f}', file=outfile)
            print('', file=outfile)

        # duplication rates
        print('# duplication rates #', file=outfile)
        print('', file=outfile)
        
        for sm in md['seq_method'].unique():
            print(f'Duplication rate (%), {sm} samples:', file=outfile)
            print(f'    Mean: {md[md["seq_method"]==sm]["pct_duplicated"].mean():.2f}', file=outfile)
            print(f'    Median: {md[md["seq_method"]==sm]["pct_duplicated"].median():.2f}', file=outfile)
            print(f'    SD: {md[md["seq_method"]==sm]["pct_duplicated"].std():.2f}', file=outfile)
            print('', file=outfile)

        # sex QC
        print('', file=outfile)
        print('## SEX QC ##', file=outfile)
        print('', file=outfile)
        
        print(f'Total number of sexed samples: {n_sexed}', file=outfile)
        print(f'    Pass sex QC: {n_pass} ({n_pass/n_sexed*100:.2f}%)', file=outfile)
        print(f'    Fail sex QC: {n_fail} ({n_fail/n_sexed*100:.2f}%)', file=outfile)
        print(f'    Suspect sex QC: {n_susp} ({n_susp/n_sexed*100:.2f}%)', file=outfile)
        print('', file=outfile)

        print(f'Total number of males: {len(males)}', file=outfile)
        print(f'    Males passing sex QC:  {len(pass_males)} ({len(pass_males)/len(males)*100:.2f}%)', file=outfile)
        print(f'    Males failing sex QC:  {len(fail_males)} ({len(fail_males)/len(males)*100:.2f}%)', file=outfile)
        print(f'    Suspect sex QC males:  {len(susp_males)} ({len(susp_males)/len(males)*100:.2f}%)', file=outfile)
        print('', file=outfile)

        print(f'Total number of females: {len(females)}', file=outfile)
        print(f'    Females passing sex QC:  {len(pass_females)} ({len(pass_females)/len(females)*100:.2f}%)', file=outfile)
        print(f'    Females failing sex QC:  {len(fail_females)} ({len(fail_females)/len(females)*100:.2f}%)', file=outfile)
        print(f'    Suspect sex QC females:  {len(susp_females)} ({len(susp_females)/len(females)*100:.2f}%)', file=outfile)
        print('', file=outfile)
        
        # genotype counts
        print('', file=outfile)
        print('## GENOTYPE COUNTS ##', file=outfile)
        print('', file=outfile)    

        # INFO scores
        print('', file=outfile)
        print('## STITCH INFO SCORES ##', file=outfile)
        print('', file=outfile)    
        print(f'Total number of SNPs genotyped: {len(info_all)}', file=outfile)
        print(f'    Mean INFO score: {info_all["INFO_SCORE"].mean():.3f}', file=outfile)
        print(f'    Median INFO score: {info_all["INFO_SCORE"].median():.3f}', file=outfile)
        print(f'    SD: {info_all["INFO_SCORE"].std():.3f}', file=outfile)
        print('', file=outfile)    
        print(f'Total number of SNPs retained (with INFO > 0.9): {len(info_clean)}', file=outfile)
        print(f'    Mean INFO score: {info_clean["INFO_SCORE"].mean():.3f}', file=outfile)
        print(f'    Median INFO score: {info_clean["INFO_SCORE"].median():.3f}', file=outfile)
        print(f'    SD: {info_clean["INFO_SCORE"].std():.3f}', file=outfile)
        print('', file=outfile)
        
        # heterozygosity
        print('', file=outfile)
        print('## HETEROZYGOSITY ##', file=outfile)
        print('', file=outfile)
        print(f'QC heterozygosity rate threshold: +/- 4 standard deviations', file=outfile)
        print(f'    Pass: het within mean +/- 4sd & SNP missingness <= {missing_rate_threshold}', file=outfile)
        print(f'    Fail: het NOT within mean +/- 4sd & SNP missingness > {missing_rate_threshold}', file=outfile)
        print(f'    Suspect: het NOT within mean +/- 4sd & SNP missingness <= {missing_rate_threshold}', file=outfile)
        print('', file=outfile)
        print(f'Total number of samples: {len(md)}', file=outfile)
        print(f'    Passed: {len(md[md["qc_heterozygosity"]=="pass"])} ({len(md[md["qc_heterozygosity"]=="pass"])/len(md)*100:.2f}%)', file=outfile)
        print(f'    Failed: {len(md[md["qc_heterozygosity"]=="fail"])} ({len(md[md["qc_heterozygosity"]=="fail"])/len(md)*100:.2f}%)', file=outfile)
        print(f'    Suspect: {len(md[md["qc_heterozygosity"]=="suspect"])} ({len(md[md["qc_heterozygosity"]=="suspect"])/len(md)*100:.2f}%)', file=outfile)
        print('', file=outfile)

        for sm in md['seq_method'].unique():
            smd = md[md['seq_method']==sm]
            print(f'Number of {sm} samples: {len(smd)}', file=outfile)
            print(f'    Passed: {len(smd[smd["qc_heterozygosity"]=="pass"])} ({len(smd[smd["qc_heterozygosity"]=="pass"])/len(smd)*100:.2f}%)', file=outfile)
            print(f'    Failed: {len(smd[smd["qc_heterozygosity"]=="fail"])} ({len(smd[smd["qc_heterozygosity"]=="fail"])/len(smd)*100:.2f}%)', file=outfile)
            print(f'    Suspect: {len(smd[smd["qc_heterozygosity"]=="suspect"])} ({len(smd[smd["qc_heterozygosity"]=="suspect"])/len(smd)*100:.2f}%)', file=outfile)
            print(f'    Mean heterozygosity rate (%): {smd["heterozygosity"].mean()*100:.2f}', file=outfile)  
            print(f'    Median heterozygosity rate (%): {smd["heterozygosity"].median()*100:.2f}', file=outfile)  
            print(f'    SD: {smd["heterozygosity"].std()*100:.2f}', file=outfile)          
            print('', file=outfile)

        print('', file=outfile)
        print('## ANALYSIS FILTERS ##', file=outfile)
        print('', file=outfile)    
        print(f'Total number of SNPs genotyped: {len(info_all)}', file=outfile)
        print(f'Total number of SNPs retained (with INFO > 0.9): {len(info_clean)}', file=outfile)
        print(f'Total number of retained SNPs with MAF > 0.01: {len(snps_maf)}', file=outfile)
        print(f'Total number of retained SNPs with missingness < 10%: {len(snps_miss)}', file=outfile)
        print(f'Total number of retained SNPs with HWE -log(P) > 6: {len(snps_hwe)}', file=outfile)
        print(f'Total number of retained SNPs with following MAF, missingness, & HWE filters: {len(snps_filtered)}', file=outfile)
        print('', file=outfile)
