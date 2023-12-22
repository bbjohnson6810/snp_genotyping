#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# Usage: python3 reads_after_demux.py -o output_file_prefix -i input_file

def help():
	print("====== plot_reads_mapped.py =====")
	print("Plot number of reads per sample after demultiplexing by Fgbio")
	print("-o <output file prefix>       the output file prefix")
	print("-i <input file>               the input file ('mapping_stats' file with samtools and mosdepth outputs)")
	print("Usage: python3 reads_after_demux.py -i input_file -o output_file_prefix")
	sys.exit()


# read in the 'mapping_stats' file and add new columns
def get_map_stats(file):
    mapped_reads = pd.read_csv(file, delimiter=',',header=0)
    mapped_reads['pct_primary_unmapped'] = 1 - mapped_reads['pct_primary_mapped']/100
    mapped_reads['duplication_rate'] = mapped_reads['duplicate_primary_reads'] / mapped_reads['primary_reads']
    return mapped_reads

# boxplots of mapped read counts, grouped by library
def plot_mapped_counts(mapped_reads, column, title, output_file):

	if len(mapped_reads['library_id'].unique()) > 4:
		plt.figure(figsize=(len(mapped_reads['library_id'].unique())*2, 8))
	else:
		plt.figure(figsize=(8, 8))

	if mapped_reads[column].mean() >= 1e6:
		mapped_reads['response'] = mapped_reads[column] / 1e6
		ax = sns.boxplot(data=mapped_reads, x='library_id', y='response')
		ax = sns.swarmplot(data=mapped_reads, x='library_id', y='response', color='.25', size=4)
		ax.set(xlabel='Library ID', ylabel='Number of Reads (millions)')
	elif 1e3 <= mapped_reads[column].mean() < 1e6:
		mapped_reads['response'] = mapped_reads[column] / 1e3
		ax = sns.boxplot(data=mapped_reads, x='library_id', y='response')
		ax = sns.swarmplot(data=mapped_reads, x='library_id', y='response', color='.25', size=4)
		ax.set(xlabel='Library ID', ylabel='Number of Reads (thousands)')
	else:
		mapped_reads['response'] = mapped_reads[column]
		ax = sns.boxplot(data=mapped_reads, x='library_id', y='response')
		ax = sns.swarmplot(data=mapped_reads, x='library_id', y='response', color='.25', size=4)
		ax.set(xlabel='Library ID', ylabel='Number of Reads')
	
	ax.set(title = title)
	plt.savefig(output_file)
	plt.close()

def plot_mapped_pct(mapped_reads, column, title, output_file):

	if len(mapped_reads['library_id'].unique()) > 4:
		plt.figure(figsize=(len(mapped_reads['library_id'].unique())*2, 8))
	else:
		plt.figure(figsize=(8, 8))

	mapped_reads['response'] = mapped_reads[column]
	ax = sns.boxplot(data=mapped_reads, x='library_id', y='response')
	ax = sns.swarmplot(data=mapped_reads, x='library_id', y='response', color='.25', size=4)
	ax.set(xlabel='Library ID', ylabel=f'Mapping Rate (%)')
	ax.set(title = title)
	plt.savefig(output_file)
	plt.close()

# parse options
if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
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
		raise "Please provide a input file ('mapped_reads' file summarizing samtools flagstat results)"


### make figures ###

# plot total reads per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'total_reads', \
	"Total Number of Reads", output_file_prefix + "_total_reads_count.png")

# plot primary reads per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'primary_reads', \
	"Number of Primary Reads", output_file_prefix + "_primary_reads_count.png")

# plot duplicate reads per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'duplicate_primary_reads', \
	"Number of Duplicate Reads", output_file_prefix + "_duplicate_reads_count.png")

# plot duplication rate per sample, grouped by library
plot_mapped_pct(get_map_stats(input_file), 'duplication_rate', \
	"Duplication Rate for Primary Reads (duplicate primary reads / all primary reads)", \
		output_file_prefix + "_mapped_reads_pct.png")

# plot total mapped reads per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'total_mapped', \
	"Total Number of Mapped Reads", output_file_prefix + "_mapped_reads_count.png")

# plot mapping rate per sample, grouped by library
plot_mapped_pct(get_map_stats(input_file), 'pct_total_mapped', \
	"Mapping Rate for All Reads (mapped reads / total reads)", \
		output_file_prefix + "_mapped_reads_pct.png")

# plot mapped primary reads per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'primary_mapped', \
	"Number of Mapped Primary Reads", output_file_prefix + "_mapped_primary_reads_count.png")

# plot mapping rate of primary reads per sample, grouped by library
plot_mapped_pct(get_map_stats(input_file), 'pct_primary_mapped', \
	"Mapping Rate for Primary Reads (mapped primary reads / all primary reads)", \
		output_file_prefix + "_mapped_primary_reads_pct.png")

# plot mapping rate of primary reads per sample, grouped by library
plot_mapped_pct(get_map_stats(input_file), 'pct_primary_unmapped', \
	"Mapping Rate for Primary Reads (unmapped primary reads / all primary reads)", \
		output_file_prefix + "_unmapped_primary_reads_pct.png")

# plot mean read depth per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'mean_read_depth', \
	"Mean Read Depth", output_file_prefix + "_read_depth_mean.png")

# plot max read depth per sample, grouped by library
plot_mapped_counts(get_map_stats(input_file), 'max_read_depth', \
	"Maximum Read Depth", output_file_prefix + "_read_depth_max.png")

