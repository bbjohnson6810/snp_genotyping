#!/usr/bin/env python3
import os
import sys
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from optparse import OptionParser
from matplotlib.ticker import NullFormatter

# Usage: python3 genotypes_pca.py -o output_file_prefix -c plink_eigenvector -v plink_eigenvalue -m metadata

def help():
	print("====== genotypes_pca.py =====")
	print("Plot Sample PCA after genotyping")
	print("-o <output file prefix>                      the output file prefix")
	print("-c <plink2 eigenvector file>            the plink2 eigenvector file")
	print("-v <plink2 eigenvalue file>              the plink2 eigenvalue file")
	print("-m <metadata file>                                the metadata file")
	print("Usage: python3 genotypes_pca.py -o output_file_prefix -c plink_eigenvector -v plink_eigenvalue -m metadata")
	sys.exit()

def read_eigenvector(file):
	eigenvec = pd.read_csv(file, delimiter="\t", dtype=str, usecols=["IID"] + ["PC"+str(i+1) for i in range(10)])
	for i in range(10):
		eigenvec["PC"+str(i+1)] = pd.to_numeric(eigenvec["PC"+str(i+1)])
	eigenvec = eigenvec.rename(columns={"IID": "sample_id"})
	return eigenvec

def read_metadata(file):
	metadata = pd.read_csv(file, delimiter=",", dtype=str, usecols=["sample_id", "rfid", "library_name", "project_name", "sex"])
	return metadata

def read_eigenvalue(file):
	eigenval_ls = pd.read_csv(file, delimiter="\t", dtype=float, header=None)[0].tolist()
	eigenval=dict()
	for i in range(10):
		eigenval["PC"+str(i+1)] = eigenval_ls[i]
	return eigenval

def plot_pca_1(eigenvec, eigenval, pc1, pc_n, hue, output_file):
	plt.figure(figsize=(8, 8))
	ax = sns.scatterplot(data=eigenvec, x=pc1, y=pc_n, hue=hue, alpha=0.5)
	ax.set(xlabel=pc1 + " (" + "{0:.2f}%".format(eigenval[pc1]/sum(eigenval.values()) * 100) + ")",
			ylabel=pc_n + " (" + "{0:.2f}%".format(eigenval[pc_n]/sum(eigenval.values()) * 100) + ")",
			title=pc1 + " vs " + pc_n + " Based on Sample Genotypes (colored by " + hue + ")")
	if len(eigenvec[hue].unique()) > 10:
		ax.legend().set_visible(False)
	plt.savefig(output_file)
	plt.close()

def plot_pca_2(eigenvec, eigenval, pc2, pc_n, hue, output_file):
	plt.figure(figsize=(8, 8))
	ax = sns.scatterplot(data=eigenvec, x=pc2, y=pc_n, hue=hue, alpha=0.5)
	ax.set(xlabel = pc2 + " (" + "{0:.2f}%".format(eigenval[pc2]/sum(eigenval.values()) * 100) + ")",
			ylabel = pc_n + " (" + "{0:.2f}%".format(eigenval[pc_n]/sum(eigenval.values()) * 100) + ")",
			title = pc2 + " vs " + pc_n + " Based on Sample Genotypes (colored by " + hue + ")")
	if len(eigenvec[hue].unique()) > 10:
		ax.legend().set_visible(False)
	plt.savefig(output_file)
	plt.close()

def plot_pca_3(eigenvec, eigenval, pc3, pc4, hue, output_file):
	plt.figure(figsize=(8, 8))
	ax = sns.scatterplot(data=eigenvec, x=pc3, y=pc4, hue=hue, alpha=0.5)
	ax.set(xlabel=pc3 + " (" + "{0:.2f}%".format(eigenval[pc3]/sum(eigenval.values()) * 100) + ")",
			ylabel=pc4 + " (" + "{0:.2f}%".format(eigenval[pc4]/sum(eigenval.values()) * 100) + ")",
			title=pc3 + " vs " + pc4 + " Based on Sample Genotypes (colored by " + hue + ")")
	if len(eigenvec[hue].unique()) > 10:
		ax.legend().set_visible(False)
	plt.savefig(output_file)
	plt.close()

# function to return how many standard deviations from the mean a given PC value is
def outlier_degree(pc_val, pc_mean, pc_std):
    pc_deg = abs(pc_val - pc_mean) / pc_std
    return pc_deg

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="out_file_prefix", help="<output file prefix>")
	parser.add_option('-c', type="string", nargs=1, dest="eigenvector", help="<plink2 eigenvector file>")
	parser.add_option('-v', type="string", nargs=1, dest="eigenvalue", help="<plink2 eigenvalue file>")
	parser.add_option('-m', type="string", nargs=1, dest="metadata", help="<metadata>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")
	options, args = parser.parse_args()
	if len(sys.argv) == 1 or options.help != None:
		help()
	if options.out_file_prefix != None:
		output_file_prefix = options.out_file_prefix
	else:
		raise "Please provide an output file"
	if options.eigenvector != None:
		eigenvector_f = options.eigenvector
	else:
		raise "Please provide a plink2 eigenvector file"
	if options.eigenvalue != None:
		eigenvalue_f = options.eigenvalue
	else:
		raise "Please provide a plink2 eigenvalue file"
	if options.metadata != None:
		metadata_f = options.metadata
	else:
		raise "Please provide a metadata file"

	eigenvector = read_eigenvector(eigenvector_f)
	metadata = read_metadata(metadata_f)
	eigenvalue = read_eigenvalue(eigenvalue_f)

	eigenvec_metadata = pd.merge(eigenvector, metadata, on=["sample_id"], how="left")

	# identify PCA outliers based on PC mean, standard deviation
	pc1_mean = eigenvec_metadata['PC1'].mean()
	pc1_std = eigenvec_metadata['PC1'].std()
	pc2_mean = eigenvec_metadata['PC2'].mean()
	pc2_std = eigenvec_metadata['PC2'].std()
	pc3_mean = eigenvec_metadata['PC3'].mean()
	pc3_std = eigenvec_metadata['PC3'].std()
	pc4_mean = eigenvec_metadata['PC4'].mean()
	pc4_std = eigenvec_metadata['PC4'].std()

	pca_outliers = eigenvec_metadata[['sample_id','library_name','project_name','sex','PC1','PC2','PC3','PC4']]
	pca_outliers['PC1_zscore'] = pca_outliers['PC1'].apply(lambda x: outlier_degree(x, pc1_mean, pc1_std))
	pca_outliers['PC2_zscore'] = pca_outliers['PC2'].apply(lambda x: outlier_degree(x, pc2_mean, pc2_std))
	pca_outliers['PC3_zscore'] = pca_outliers['PC3'].apply(lambda x: outlier_degree(x, pc3_mean, pc3_std))
	pca_outliers['PC4_zscore'] = pca_outliers['PC4'].apply(lambda x: outlier_degree(x, pc4_mean, pc4_std))

	# save csv of PCA values and standard deviations
	pca_outliers.to_csv(output_file_prefix + '_pca_vals.csv', index=False)

	extreme_outliers = pca_outliers[(pca_outliers['PC1_zscore'] > 5) | (pca_outliers['PC2_zscore'] > 5) | \
								(pca_outliers['PC3_zscore'] > 5) | (pca_outliers['PC4_zscore'] > 5)]

	pca_outliers.to_csv(output_file_prefix + '_pca_vals.csv', index=False)
	extreme_outliers.to_csv(output_file_prefix + '_pca_outliers.csv', index=False)
	
	# PCA plots
	for i in [2,3,4]:
		for hue in ["library_name", "project_name", "sex"]:
			plot_pca_1(eigenvec_metadata, eigenvalue,
					"PC1", f"PC{str(i)}", hue, 
					output_file_prefix + "_sample_PC1_vs_PC" + str(i) + "_" + hue + ".png")

	for i in [3,4]:
		for hue in ["library_name", "project_name", "sex"]:
			plot_pca_2(eigenvec_metadata, eigenvalue,
					"PC2", f"PC{str(i)}", hue, 
					output_file_prefix + "_sample_PC2_vs_PC" + str(i) + "_" + hue + ".png")

	for hue in ["library_name", "project_name", "sex"]:
		plot_pca_3(eigenvec_metadata, eigenvalue,
				"PC3", "PC4", hue, 
				output_file_prefix + "_sample_PC3_vs_PC4_" + hue + ".png")
