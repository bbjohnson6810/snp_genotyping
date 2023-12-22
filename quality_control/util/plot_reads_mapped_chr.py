#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 reads_after_mkDup_chr.py -s metadata_file -i input_file -o output_file_prefix

def help():
	print("====== reads_after_mkDup_chr.py =====")
	print("Plot number of mapped reads per sample per chromosome after alignment by bwa and mark duplicates by picard")
	print("-s <metadata file>                                            the metadata file")
	print("-o <output file prefix>                                  the output file prefix")
	print("-i <input file>                      the input file (mapped reads per chr file)")
	print("Usage: python3 reads_after_mkDup_chr.py -i input_file -o output_file_prefix")
	sys.exit()

def read_sample_sheet(metadata_file):
	sample_sheet = pd.read_csv(metadata_file, delimiter=",", dtype=str)
	sample_sheet = sample_sheet.sort_values(by=['sample_id']).reset_index(drop=True)
	return sample_sheet


def read_mapped_chr(file):
	mapped_chr = pd.read_csv(file, delimiter=',', header=0)
	mapped_chr = mapped_chr.sort_values(by=['sample_id']).reset_index(drop=True)
	chrs = list(range(1,21)) + ['X', 'Y', 'M']
	chr_cols = ['chr' + str(chr) + '_mapped' for chr in chrs]
	pct_chrs = []  
	for col in chr_cols:
		chr = col.split('_')[0]
		mapped_chr[chr + '_pct_mapped'] = mapped_chr[col] / mapped_chr['total_reads']
		pct_chrs.append(chr + '_pct_mapped')
	mapped_chr = mapped_chr.sort_values(by=['sample_id']).reset_index(drop=True)
	return mapped_chr, chr_cols, pct_chrs


def plot_mapped_reads_chr(mapped_chr, chr_cols, output_file):
	
	chr_list = list(range(1,21)) + ['X', 'Y', 'M']
	chr_col_list = list(mapped_chr[chr_cols].columns)
	mapped_chr.rename(columns = dict(zip(chr_col_list,chr_list)), inplace=True)
	chr_cols = chr_list

	chr_df = pd.melt(mapped_chr[["sample_id", "library_id"] + chr_cols], 
			 id_vars=["sample_id", "library_id"], value_vars=chr_cols,var_name="chr", value_name="mapped_reads")
	chr_df['mapped_reads'] = chr_df['mapped_reads'] / 1e6

	plt.figure(figsize=(15, 6))
	if len(mapped_chr["library_id"].unique()) > 10:
		ax = sns.boxplot(data = chr_df, x="chr", y="mapped_reads")
		ax.set(xlabel="Chromosome", ylabel="Mapped Reads (million)",
		title="Number of Mapped Reads per Chromosome (across all libraries)")
	else:
		ax = sns.boxplot(data = chr_df, x="chr", y="mapped_reads", hue="Library_ID")
		ax.set(xlabel="Chromosome", ylabel="Mapped Reads (million)",
		title="Number of Mapped Reads per Chromosome (across all libraries)")
	plt.savefig(output_file)
	plt.close()
	return chr_df

def plot_percent_mapped_chr(mapped_chr, pct_chrs, output_file):
	
	chr_list = list(range(1,21)) + ['X', 'Y', 'M']
	chr_col_list = list(mapped_chr[pct_chrs].columns)
	mapped_chr.rename(columns = dict(zip(chr_col_list,chr_list)), inplace=True)
	chr_cols = chr_list

	chr_df = pd.melt(mapped_chr[["sample_id", "library_id"] + chr_cols], 
			 id_vars=["sample_id", "library_id"], value_vars=chr_cols,var_name="chr", value_name="percent_mapped")
	
	plt.figure(figsize=(15, 6))
	if len(mapped_chr["library_id"].unique()) > 10:
		ax = sns.boxplot(data=chr_df, x="chr", y="percent_mapped")
	else:
		ax = sns.boxplot(data=chr_df, x="chr", y="percent_mapped", hue="library_id")
	vals = ax.get_yticks()
	ax.set(xlabel="Chromosome", ylabel="Percentage Mapped Reads",
		title="Percentage Mapped Reads per Chromosome")
	plt.savefig(output_file)
	plt.close()


if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-s', type="string", nargs=1, dest="metadata_file", help="<metadata file>")
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-i', type="string", nargs=1, dest="in_file", help="<input file>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.in_file != None:
		input_file = options.in_file
	else:
		raise "Please provide a input file (mapped reads per chr file)"
	if options.metadata_file != None:
		metadata_file = options.metadata_file
	else:
		raise "Please provide a metadata file with Sample_ID and sex info"

	mapped_chr, chr_cols, pct_chrs = read_mapped_chr(input_file)
	sample_sheet = read_sample_sheet(metadata_file)


	# create figures
	plot_mapped_reads_chr(mapped_chr, chr_cols, output_file_prefix + "_mapped_reads.png")
	plot_percent_mapped_chr(mapped_chr, pct_chrs, output_file_prefix + "_percent_mapped.png")
