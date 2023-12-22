#!/usr/bin/env python3

import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
import seaborn as sns
from optparse import OptionParser

def help():
	print("====== plot_sample_het_vs_missing.py =====")
	print("Plot Sample heterozygosity vs. missingness")
	print("-i <info_scores>            file of STITCH INFO scores")
	print("-a <afreq>                  plink allele frequencies file")
	print("-w <hwe>                    plink HWE file")
	print("-x <vmiss>                  plink genotype missingness file")
	print("-o <output file prefix>     the output file prefix")
	print("Usage: python3 genotypes_pca.py -i gtype_log -o out_prefix")
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

# function to plot snp stats
# note that plotted threshold lines are hypothetical: we don't actually apply a missingness or HWE threshold
# but these plots show what would be included in the dataset if we did apply these cutoffs
def plot_snps(snp_info, x, y, xlabel, ylabel, 
              MAF_cutoff=None, hwe_cutoff=None, missing_cutoff=None, outfile=''):

    nullfmt = NullFormatter()
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_hist = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]

    # start with a rectangular figure
    plt.figure(1, figsize=(8.5, 8.5), dpi=130)
    axHist = plt.axes(rect_hist)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)

    # the main plot:
    sns.histplot(ax=axHist, data=snp_info, bins=50, x=x, y=y)
    axHist.set(xlabel=xlabel, ylabel=ylabel)
    if hwe_cutoff is not None:
        axHist.axhline(y=hwe_cutoff, color='red', linestyle='--', label=f'Hardy-Weinberg -log(P) threshold: {hwe_cutoff}')
    if missing_cutoff is not None:
        axHist.axhline(y=missing_cutoff, color='red', linestyle='--', label=f'SNP missingness threshold: {missing_cutoff}')
    if MAF_cutoff is not None:
        axHist.axvline(MAF_cutoff,color='orange',ls='--')

    # top plot
    sns.histplot(ax=axHistx, data=snp_info, x=x, bins=50)
    if MAF_cutoff is not None:
        axHistx.axes.set_title(f'{ylabel} vs MAF > {MAF_cutoff:.3f} ({len(snp_info)} SNPs)',fontsize=16)
        axHistx.axvline(MAF_cutoff,color='orange',ls='--')

    else:
        axHistx.axes.set_title(f'{ylabel} vs {xlabel} ({len(snp_info)} SNPs)',fontsize=16)
    axHistx.set_xlabel('',fontsize=20)
    axHistx.set_ylabel('Number of SNPs')
    axHistx.tick_params(labelsize=5)

    # right plot
    sns.histplot(ax=axHisty, data=snp_info, y=y, bins=50)
    axHisty.set(xlabel='Number of SNPs', ylabel='',title='')
    if hwe_cutoff is not None:
        axHisty.axhline(hwe_cutoff,color='red',ls='--')
    elif missing_cutoff is not None:
        axHisty.axhline(missing_cutoff,color='red',ls='--')
    
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


# options parser
if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-i', type="string", nargs=1, dest="info_file", help="<info scores file>")
    parser.add_option('-a', type="string", nargs=1, dest="afreq_file", help="<plink .afreq file>")
    parser.add_option('-w', type="string", nargs=1, dest="hwe_file", help="<plink .hardy file>")
    parser.add_option('-x', type="string", nargs=1, dest="vmiss_file", help="<plink .vmiss file")
    parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()
    
    if len(sys.argv) == 1 or options.help != None:
        help()
    if options.out_file_prefix != None:
        output_file_prefix = options.out_file_prefix
    else:
        raise "Please provide an output prefix"
    if options.info_file != None:
        info = options.info_file
    else:
        raise "Please provide a file of STITCH INFO scores"
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


    # create snp dataframe
    snp_df = make_snp_df(info, afreq, hwe, vmiss)

    print('Plotting...')

    # plot HWE ~ MAF for all snps
    plot_snps(snp_df, 'MAF', '-log(P)', 'Minor Allele Frequency', 'Hardy-Weinberg -log(P)', hwe_cutoff=10,
            outfile = output_file_prefix + '_hwe_vs_maf.png')
    
    # plot HWE ~ MAF for snps with MAF > 0
    plot_snps(snp_df[snp_df['MAF'] > 0].reset_index(drop=True),
                'MAF', '-log(P)', 'Minor Allele Frequency', 'Hardy-Weinberg -log(P)', MAF_cutoff=0, hwe_cutoff=10,
            outfile = output_file_prefix + '_hwe_vs_maf0.png')
    
    # plot HWE ~ MAF for snps with MAF > 0.005
    plot_snps(snp_df[snp_df['MAF'] > 0.005].reset_index(drop=True),
                'MAF', '-log(P)', 'Minor Allele Frequency', 'Hardy-Weinberg -log(P)', MAF_cutoff=0.005, hwe_cutoff=10,
            outfile = output_file_prefix + '_hwe_vs_maf0.005.png')
    
    # plot snp missing rate ~ MAF for all snps
    plot_snps(snp_df, 'MAF', 'missing_rate', 'Minor Allele Frequency', 'Missing Rate', missing_cutoff=0.1,
            outfile = output_file_prefix + '_missing_vs_maf.png')
    
    # plot snp missing rate ~ MAF for snps with MAF > 0
    plot_snps(snp_df[snp_df['MAF'] > 0].reset_index(drop=True),
                'MAF', 'missing_rate', 'Minor Allele Frequency', 'Missing Rate', MAF_cutoff=0, missing_cutoff=0.1,
                outfile = output_file_prefix + '_missing_vs_maf0.png')
    
    # plot snp missing rate ~ MAF for snps with MAF > 0.005
    plot_snps(snp_df[snp_df['MAF'] > 0.005].reset_index(drop=True),
                'MAF', 'missing_rate', 'Minor Allele Frequency', 'Missing Rate', MAF_cutoff=0.005, missing_cutoff=0.1,
            outfile = output_file_prefix + '_missing_vs_maf0.005.png')