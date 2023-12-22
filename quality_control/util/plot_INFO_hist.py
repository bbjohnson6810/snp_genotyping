#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser

# use this script to plot a histogram of per-snp INFO scores produced by STITCH
# Usage: python3 genotypes_imputation_INFO.py -o output_file_prefix -i input_file

def help():
	print("====== genotypes_imputation_INFO.py =====")
	print("Plot SNPs INFO score after genotypes imputation")
	print("-o <output file prefix>                              the output file prefix")
	print("-i <input file>                                 the input file (INFO score)")
	print("Usage: python3 genotypes_imputation_INFO.py -i input_file -o output_file_prefix")
	sys.exit()

# function to read in column of INFO scores
def read_STITCH_INFO(file):
	STITCH_INFO = pd.read_csv(file, delimiter="\t", dtype=str)
	STITCH_INFO["INFO_SCORE"] = pd.to_numeric(STITCH_INFO["INFO_SCORE"])
	return STITCH_INFO

# function to plot a histogram of INFO scores
def plot_INFO_hist(STITCH_INFO, INFO_threshold, output_file):
	plt.figure(figsize=(8, 8))
	ax = sns.histplot(data=STITCH_INFO, x="INFO_SCORE", kde=True, bins=80)
	ax.axvline(x=INFO_threshold, color="red", linestyle="--",
		label="INFO score threshold: " +str(INFO_threshold)+ " (" +
		str(len(STITCH_INFO[STITCH_INFO["INFO_SCORE"] > INFO_threshold])) +" SNPs)")
	ax.set(xlabel="INFO Score", ylabel="Count", title=f"SNPs INFO Score ({len(STITCH_INFO)} snps)")
	ax.legend()
	plt.savefig(output_file)
	plt.close()

# options parser
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
		raise "Please provide a input file (INFO score)"

# plot a histogram of INFO scores passing a given threshold score
	INFO_threshold = 0.9
	plot_INFO_hist(read_STITCH_INFO(input_file), INFO_threshold, output_file_prefix + "_hist.png")

