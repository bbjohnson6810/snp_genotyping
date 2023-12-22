#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# use this script to 
# Usage: python3 genotypes_missing_vs_mapped_reads.py -o output_file_prefix -r mkDup_metrics -m plink_sample_missing

# help function
def help():
	print("====== genotypes_missing_vs_mapped_reads.py =====")
	print("Plot sample missing rate vs mapped reads after genotyping")
	print("-o <output file prefix>                                     the output file prefix")
	print("-r <Picard mkDup metrics file>                       the Picard mkDup metrics file")
	print("-m <plink2 sample-based missing data>         the plink2 sample-based missing data")
	print("Usage: python3 genotypes_missing_vs_mapped_reads.py -o output_file_prefix -r mkDup_metrics -m plink_sample_missing")
	sys.exit()

# function to read a file of mapping statistics
def read_mkDup_metrics(file):
	mkDup_metrics = pd.read_csv(file, delimiter="\t", dtype=str,usecols=["Sample_ID", "READ_PAIRS_EXAMINED", "UNPAIRED_READS_EXAMINED"])
	mkDup_metrics["Library_ID"] = mkDup_metrics["Sample_ID"].apply(lambda x: '_'.join(x.split('_')[:-1]))
	mkDup_metrics["READ_PAIRS_EXAMINED"] = pd.to_numeric(mkDup_metrics["READ_PAIRS_EXAMINED"])
	mkDup_metrics["UNPAIRED_READS_EXAMINED"] = pd.to_numeric(mkDup_metrics["UNPAIRED_READS_EXAMINED"])
	mkDup_metrics["MAPPED_READS"] = mkDup_metrics["UNPAIRED_READS_EXAMINED"] + mkDup_metrics["READ_PAIRS_EXAMINED"]*2
	mkDup_metrics["MAPPED_READS"] = mkDup_metrics["MAPPED_READS"]/1e6
	return mkDup_metrics

def read_mapstats(file):
	mapstats = pd.read_csv(file, usecols=['sample_id', 'library_id', 'primary_reads', 'primary_mapped'])
	mapstats['primary_reads'] = pd.to_numeric(mapstats['primary_reads'])
	mapstats['primary_mapped'] = pd.to_numeric(mapstats['primary_mapped'])
	mapstats['primary_reads'] = mapstats['primary_reads'] / 1e6
	mapstats['primary_mapped'] = mapstats['primary_mapped'] / 1e6
	return mapstats

def read_sample_missing(file):
	sample_missing = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["IID", "F_MISS"])
	sample_missing["F_MISS"] = pd.to_numeric(sample_missing["F_MISS"])
	sample_missing = sample_missing.rename(columns={"IID": "sample_id", "F_MISS": "sample_missing_rate"})
	return sample_missing

def plot_sample_missing_vs_mapped(mapstats_file, sample_missing_rate_threshold, output_file):
	nullfmt = NullFormatter()
	# definitions for the axes
	left, width = 0.1, 0.65
	bottom, height = 0.1, 0.65
	bottom_h = left_h = left + width + 0.02
	rect_hist = [left, bottom, width, height]
	rect_histx = [left, bottom_h, width, 0.2]
	rect_histy = [left_h, bottom, 0.2, height]
	# start with a rectangular Figure
	plt.figure(1, figsize=(10, 10))
	axHist = plt.axes(rect_hist)
	axHistx = plt.axes(rect_histx)
	axHisty = plt.axes(rect_histy)
	axHistx.xaxis.set_major_formatter(nullfmt)
	axHisty.yaxis.set_major_formatter(nullfmt)
	# the main plot:
	sns.scatterplot(ax=axHist, data=mapstats_file,
						x="primary_mapped", y="sample_missing_rate", alpha=0.5)
	axHist.set(xlabel="Mapped Reads (million)", ylabel="Missing Rate")
	axHist.axhline(y=sample_missing_rate_threshold, color="red", linestyle="--", label="Missing Rate Threshold: " +
		str(sample_missing_rate_threshold) + " (" + str(len(sample_missing_mkDup_metrics[sample_missing_mkDup_metrics["QC_sample_missing_rate"] == "fail"]))+
		" samples)")
	axHist.legend()
	# sub plots
	sns.histplot(ax=axHistx, data=mapstats_file,
				x="primary_mapped", bins=100, kde=True)
	axHistx.set(xlabel="", ylabel="Number of Samples", title="Sample Missing Rate vs Mapped Reads")
	sns.histplot(ax=axHisty, data=sample_missing_mkDup_metrics, y="sample_missing_rate", bins=100, kde=True)
	axHisty.set(xlabel="Number of Samples", ylabel="",title="")
	plt.savefig(output_file)
	plt.close()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-r', type="string", nargs=1, dest="mapstats_file", help="<Picard mkDup metrics file>")
	parser.add_option('-m', type="string", nargs=1, dest="smiss", help="<plink2 sample-based missing data")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide a output file"
	if options.picard_mkDup_metrics != None:
		picard_mkDup_metrics = options.picard_mkDup_metrics
	else:
		raise "Please provide a mapstats file"
	if options.smiss != None:
		smiss = options.smiss
	else:
		raise "Please provide a plink2 sample-based missing data"

	# read in mapping statistics from Picard
	mapstats = read_mapstats(mapstats_file)
	sample_missing = read_sample_missing(smiss)

	sample_missing_mkDup_metrics = pd.merge(sample_missing, mapstats, on=["sample_id"], how="left")
	sample_missing_rate_threshold = 0.1
	sample_missing_mkDup_metrics["QC_sample_missing_rate"] = sample_missing_mkDup_metrics["sample_missing_rate"].apply(lambda x: "pass" if x < sample_missing_rate_threshold else "fail")

	plot_sample_missing_vs_mapped(sample_missing_mkDup_metrics, sample_missing_rate_threshold, output_file_prefix + "_sample_missing_vs_mapped_reads.png")

	QC_sample_missing_rate_threshold_10 = sample_missing_mkDup_metrics[["sample_id", "library_id", "sample_missing_rate", "QC_sample_missing_rate"]]
	QC_sample_missing_rate_threshold_10.to_csv(output_file_prefix+"QC_sample_missing_rate_threshold_10pct.csv", sep=',', index=False)