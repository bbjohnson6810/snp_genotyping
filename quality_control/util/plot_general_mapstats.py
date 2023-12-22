#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 plot_general_mapstats.py 

def help():
    print("====== reads_after_mkDup_chr.py =====")
    print("Plot number of mapped reads per sample per chromosome after alignment by bwa and mark duplicates by picard")
    print("-s <metadata file>                                            the metadata file")
    print("-o <output file prefix>                                  the output file prefix")
    print("-i <genotype log>                      the general genotyping log for all samples")
    print("Usage: python3 reads_after_mkDup_chr.py -i input_file -o output_file_prefix")
    sys.exit()

 
def plot_mapping_by_seq_method_violin(md, output_file):

    plt.figure(figsize=(8, 6), dpi = 130)
    ax = sns.violinplot(data=md, x='primary_mapped', y='seq_method')
    ax.set(ylabel='Samples', xlabel='Number of Mapped Reads (million)',
        title='Number of Reads Mapped')
    plt.xticks(rotation=45)
    plt.savefig(output_file)
    plt.close()

def plot_mapping_by_seq_method_hist(md, output_file):

    plt.figure(figsize=(6, 5), dpi = 130)
    ax = sns.histplot(data=md,
                    x='primary_mapped', hue='seq_method', bins=100)
    ax.set(ylabel='Number of Samples', xlabel='Number of Mapped Reads (million)')
    plt.savefig(output_file)
    plt.close()


def plot_mapping_rate_by_seq_method_violin(md, output_file):
	
    # calculate mapping rates
    md['pct_mapped'] = md['primary_mapped'] / md['primary_reads'] * 100

    plt.figure(figsize=(8, 6), dpi = 130)
    ax = sns.violinplot(data=md, x='pct_mapped', y='seq_method')
    ax.set(ylabel='Samples', xlabel='Mapping Rate (%)', title='Mapping Rate')
    plt.xticks(rotation=45)
    plt.savefig(output_file)
    plt.close()

def plot_mapping_rate_by_seq_method_hist(md, output_file):
	
    # calculate mapping rates
    md['pct_mapped'] = md['primary_mapped'] / md['primary_reads'] * 100

    plt.figure(figsize=(8, 6), dpi = 130)
    ax = sns.histplot(data=md, x='pct_mapped', hue='seq_method', bins=100)
    ax.set(ylabel='Samples', xlabel='Mapping Rate (%)', title='Mapping Rate')
    plt.savefig(output_file)
    plt.close()

def plot_pct_duplicate_reads(md, output_file):
	
    # calculate mapping rates
    md['pct_duplicated'] = md['duplicate_primary_reads'] / md['primary_reads'] * 100
    
    plt.figure(figsize=(8, 6), dpi = 130)
    ax = sns.violinplot(data=md, x='pct_duplicated', y='seq_method')
    ax.set(ylabel='Samples', xlabel='Read Duplication Rate (%)', title='Percentage of Reads Duplicated')
    plt.xticks(rotation=45)
    plt.savefig(output_file)
    plt.close()


if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-i', type="string", nargs=1, dest="gtype_log", help="<genotype log>")
    parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
    parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
    options, args = parser.parse_args()

    print(f"gtype_log: {options.gtype_log}")
    print(f"out_file_prefix: {options.out_file_prefix}")

    if len(sys.argv) == 1 or options.help != None:
        help()
    if options.gtype_log != None:
        gtype_log = options.gtype_log
    else:
        raise "Please provide a genotyping log"
    if options.out_file_prefix != None:
        output_file_prefix = options.out_file_prefix
    else:
        raise Exception("Please provide an output prefix")
 
    # read in metadata
    gtype_log = pd.read_csv(gtype_log)

    # create figures
    plot_mapping_by_seq_method_violin(gtype_log, output_file_prefix + '_reads_mapped_seq_method_violin')
    plot_mapping_by_seq_method_hist(gtype_log, output_file_prefix + '_reads_mapped_seq_method_hist')
    plot_mapping_rate_by_seq_method_violin(gtype_log, output_file_prefix + '_mapping_rate_seq_method_violin')
    plot_mapping_rate_by_seq_method_hist(gtype_log, output_file_prefix + '_mapping_rate_seq_method_hist')
    plot_pct_duplicate_reads(gtype_log, output_file_prefix + '_duplicate_reads_pct')