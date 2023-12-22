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
	print("-i <genotype log>            the genotyping log")
	print("-o <output file prefix>                      the output file prefix")
	print("Usage: python3 genotypes_pca.py -i gtype_log -o out_prefix")
	sys.exit()

# function to read in genotype log
def read_geno_log(file):
      
      geno_log = pd.read_csv(file)
      return geno_log


# function to plot heterozygosity QC for ddGBS samples
def plot_sample_het_vs_missing(geno_log, seq_method, missing_rate_threshold, het_mean, het_std, upper_std, lower_std, outfile='', title=''):

        
    nullfmt = NullFormatter()
    # definitions for the axes
    left, width = 0.1, 0.65
    bottom, height = 0.1, 0.65
    bottom_h = left_h = left + width + 0.02
    rect_hist = [left, bottom, width, height]
    rect_histx = [left, bottom_h, width, 0.2]
    rect_histy = [left_h, bottom, 0.2, height]
    # start with a rectangular Figure
    plt.figure(1, figsize=(7, 7), dpi = 130)
    axHist = plt.axes(rect_hist)
    axHistx = plt.axes(rect_histx)
    axHisty = plt.axes(rect_histy)
    axHistx.xaxis.set_major_formatter(nullfmt)
    axHisty.yaxis.set_major_formatter(nullfmt)
    
    # the main plot:
    sns.scatterplot(ax=axHist, data=geno_log, x='missing_rate', y='heterozygosity', alpha=0.5)
    axHist.set(xlabel='Missing Rate', ylabel='Heterozygosity Rate')
    axHist.axhline(y=het_mean + upper_std*het_std, color='orange', linestyle='--', 
                   label= f'Heterozygosity Rate: {het_mean + upper_std*het_std:.3f} (Mean + {upper_std} STD)')
    axHist.axhline(y=het_mean - lower_std*het_std, color='orange', linestyle='--', 
                   label= f'Heterozygosity Rate: {het_mean - lower_std*het_std:.3f} (Mean - {lower_std} STD)')
    axHist.axvline(x=missing_rate_threshold, color='red', linestyle='--', label='Missing Rate: ' + str(missing_rate_threshold))
    axHist.legend()
        
    # sub plots
    sns.histplot(ax=axHistx, data=geno_log, x='missing_rate', bins=100, kde=True)
    axHistx.set(xlabel='', ylabel='# of Samples', title='Sample Heterozygosity Rate vs Missing Rate (' + seq_method + ')')
    # axHistx.set_xlim(-5, 135)
    axHistx.axvline(missing_rate_threshold, color='red', linestyle='--')
    sns.histplot(ax=axHisty, data=geno_log, y='heterozygosity', bins=100, kde=True)
    axHisty.set(xlabel='# of Samples', ylabel='',title='')
    axHisty.axhline(het_mean+upper_std*het_std, color='orange', linestyle='--')
    axHisty.axhline(het_mean-lower_std*het_std, color='orange', linestyle='--')
    plt.savefig(outfile, bbox_inches='tight')
    plt.close()


# options parser
if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-i', type="string", nargs=1, dest="gtype_log", help="<genotype log>")
    parser.add_option('-m', type="string", nargs=1, dest="seq_method", help="<sequencing method>")
    parser.add_option('-u', type="string", nargs=1, dest="upper_std_cutoff", help="<upper het std cutoff>")
    parser.add_option('-l', type="string", nargs=1, dest="lower_std_cutoff", help="<lower het std cutoff>")
    parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()
    if len(sys.argv) == 1 or options.help != None:
        help()
    if options.seq_method != None:
        seq_method = options.seq_method
    else:
        raise "Please provide a sequencing method (ddGBS or lcWGS)"
    if options.upper_std_cutoff != None:
        upper_std = options.upper_std_cutoff
    else:
        raise "Please provide an upper heterozygosity rate cutoff (in standard deviations)"
    if options.lower_std_cutoff != None:
        lower_std = options.lower_std_cutoff
    else:
        raise "Please provide a lower heterozygosity rate cutoff (in standard deviations)"
    if options.out_file_prefix != None:
        output_file_prefix = options.out_file_prefix
    else:
        raise "Please provide a output file"
    if options.gtype_log != None:
        gtype_log = options.gtype_log
    else:
        raise "Please provide a input file (genotype log)"


    # read in genotype log
    geno_log = read_geno_log(gtype_log)

    # subset by sequencing method
    geno_log = geno_log[geno_log['seq_method'] == seq_method]
    # ddGBS_log = geno_log[geno_log['seq_method']=='ddGBS']
    # lcWGS_log = geno_log[geno_log['seq_method']=='lcWGS']

    # missing rate threshold used in missing rate QC
    missing_rate_threshold = 0.1

    # heterozygosity summary stats based on sequencing method
    het_mean = geno_log['heterozygosity'].mean()
    het_std = geno_log['heterozygosity'].std()
    # het_mean_GBS = ddGBS_log['heterozygosity'].mean()
    # het_std_GBS = ddGBS_log['heterozygosity'].std()
    # het_mean_WGS = lcWGS_log['heterozygosity'].mean()
    # het_std_WGS = lcWGS_log['heterozygosity'].std()


    # save plots
    plot_sample_het_vs_missing(geno_log, seq_method , missing_rate_threshold, 
                               het_mean, het_std, outfile = output_file_prefix + '.png')
    