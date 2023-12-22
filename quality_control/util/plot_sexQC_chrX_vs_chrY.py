#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 reads_after_mkDup_chr.py -s metadata_file -i input_file -o output_file_prefix

def help():
	print("====== reads_after_mkDup_chr.py =====")
	print("Plot number of mapped reads per sample per chromosome after alignment by bwa and mark duplicates by picard")
	print("-i <genotype log>                      the genotype log for all samples")
	print("-o <output file prefix>                                  the output file prefix")
	print("Usage: python3 reads_after_mkDup_chr.py -i input_file -o output_file_prefix")
	sys.exit()


def plot_QC_sex(md, output_file):
    md = pd.read_csv(md)
    plt.figure(figsize=(8, 8))
    ax = sns.scatterplot(data=md, x='pct_chrX', y='pct_chrY', hue='sex', palette={'M': 'C0', 'F': 'C1'}, 
    					alpha=0.5, style='qc_sex', style_order=['pass', 'fail', 'suspect'])
    x = np.linspace(0,10,100)
    y = x/20.5+0.07
    plt.plot(x, y, 'black', linestyle='-', alpha=0.8)
    plt.plot(x, y+0.02, c='red', linestyle='--', alpha=0.8)
    plt.plot(x, y-0.02, c='red', linestyle='--', alpha=0.8)
    ax.axes.set_title('Percentage of  Mapped Reads on ChrY vs. ChrX',
    					fontsize=12)
    ax.set_ylabel('Percentage of Mapped Reads on ChrY (%)',fontsize=12)
    ax.set_xlabel('Percentage of Mapped Reads on ChrX (%)',fontsize=12)
    ax.tick_params(labelsize=12)
    plt.savefig(output_file, bbox_inches='tight')
    plt.close()

if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-i', type="string", nargs=1, dest="geno_log", help="<genotype log>")
    parser.add_option('-o', type="string", nargs=1, dest="out_prefix", help="<output file prefix>")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()
       
    if len(sys.argv) == 1 or options.help != None:
        help()
    if options.geno_log != None: 
        gtype_log = options.geno_log
    else:
        raise "Please provide a input file (genotype log)"
    if options.out_prefix != None: 
        output_file_prefix = options.out_prefix
    else:
        raise "Please provide an output file prefix"

    md = pd.read_csv(gtype_log)

    # plot figure
    plot_QC_sex(gtype_log, output_file_prefix + '_chrX_vs_chrY.png')
