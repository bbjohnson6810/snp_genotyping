#!/usr/bin/env python3
import os
import sys
import pandas as pd
import numpy as np
from datetime import datetime
from optparse import OptionParser

# use this script to request metadata either for a new sequencing pool 
# (a new flowcell to be included in an upcoming genotyping run)
# or to request updated metadata for samples from a previous run prior to re-running in a new genotyping round
# usage: python3 request_metadata.py -o output_prefix -m previous_metadata -g previous_genotype_log
# or:    python3 request_metadata.py -o output_prefix -n new_pool -r flowcell_id

def help():
	print("====== request_metadata.py =====")
	print("Request metadata for all samples in this genotyping run")
	print("-o <output file prefix>         the output file prefix")
	print("-m <metadata file>              the metadata file from the previous genotyping round (optional)")
	print("-g <genotyping log>             the complete genotyping log from the previous genotyping round (optional)")
	print("-n <new pool>                   a flowcell's 'sample barcode list' with sequencing metadata")
	print("-f <flowcell ID>                (str) the flowcell ID or 'run ID' used to produce the new pool of data")
	print("Usage 1: python3 request_metadata.py -o output_prefix -m previous_metadata -g previous_genotype_log")
	print("Usage 2: python3 request_metadata.py -o output_prefix -n new_pool -r flowcell_id")
	sys.exit()

if __name__=="__main__":
	usage = "usage: %prog [options]"
	parser = OptionParser(usage)
	parser.add_option('-o', type="string", nargs=1, dest="outfile_prefix", help="<output file prefix>")
	parser.add_option('-m', type="string", nargs=1, dest="metadata_file", help="<previous metadata file>")
	parser.add_option('-g', type="string", nargs=1, dest="gtype_log", help="<previous genotyping log>")
	parser.add_option('-n', type="string", nargs=1, dest="new_pool", help="<sample barcode list>")
	parser.add_option('-f', type="string", nargs=1, dest="flowcell_id", help="<flowcell run ID>")
	parser.add_option('-H', action="store_true", dest="help", help="Displays help screen")

	options, args = parser.parse_args()

	if len(sys.argv) == 1 or options.help != None:
		help()

	if options.outfile_prefix != None:
		outfile_prefix = options.outfile_prefix
	else:
		raise 'Please provide an output file prefix'

	if options.metadata_file != None and options.gtype_log == None:
		raise 'Please provide a genotyping log with the previous metadata file'
	
	if options.metadata_file == None and options.gtype_log != None:
		raise 'Please provide a previous metadata file with the genotyping log'

	#### method for previous metadata + genotype log
	if options.metadata_file != None and options.gtype_log != None:

		# read in previous metadata
		prev_metadata = pd.read_csv(options.metadata_file, dtype=str)

		# read in previous genotyping log
		prev_log = pd.read_csv(options.gtype_log, dtype=str)

		# keep passing samples from the previous log
		keep_categories = ['analysis', 'keep']
		keep_previous = prev_log[prev_log['sample_use'].isin(keep_categories)]['rfid'].tolist()
		keep_metadata = prev_metadata[prev_metadata['rfid'].isin(keep_previous)]

		# set metadata columns to request
		cols_to_request = ['rfid', 'library_name', 'project_name', 'flowcell_id', 'barcode', 'pcr_barcode', 
		'organism', 'strain', 'sex', 'coatcolor', 'sires', 'dams', 'fastq_files']

		# rename columns (if needed)
		keep_metadata.rename({'library':'library_name', 'runid':'flowcell_id'}, axis=1, inplace=True)

		# keep only desired columns (if needed)
		keep_metadata = keep_metadata[cols_to_request]

		# save to file
		keep_metadata.to_csv(f'{outfile_prefix}_request_new_metadata.csv', index=False)



	#### method for a new sample barcode list
	if options.metadata_file == None and options.gtype_log == None and options.new_pool != None:

		if options.flowcell_id == None:
			raise 'Please provide the flowcell ID for the sequencing pool'

		# read in sample barcode list
		new_pool = pd.read_csv(options.new_pool, dtype=str)

		# set metadata columns to request
		cols_to_request = ['rfid', 'library_name', 'project_name', 'flowcell_id', 'barcode', 'pcr_barcode', 
		'organism', 'strain', 'sex', 'coatcolor', 'sires', 'dams', 'fastq_files']

		# rename columns (if needed)
		new_pool.rename({'SampleID':'rfid', 'Barcode':'barcode', 'PCR.Barcode':'pcr_barcode',
			'PCR Barcode':'pcr_barcode', 'Library':'library_name', 'Flowcell Lane':'flowcell_id', 
			'Flowcell.Lane':'flowcell_id'}, axis=1, inplace=True)
		
		# add new columns
		new_pool['flowcell_id'] = options.flowcell_id
		new_pool['project_name'] = np.nan
		new_pool['organism'] = np.nan
		new_pool['strain'] = np.nan
		new_pool['sex'] = np.nan
		new_pool['coatcolor'] = np.nan
		new_pool['sires'] = np.nan
		new_pool['dams'] = np.nan
		new_pool['fastq_files'] = np.nan

		# keep only desired columns 
		new_pool = new_pool[cols_to_request]

		# save to file
		new_pool.to_csv(f'{outfile_prefix}_request_metadata.csv', index=False)
	