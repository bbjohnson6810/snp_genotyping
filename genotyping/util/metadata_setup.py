import pandas as pd
import sys

# use this script to read in a flowcell sample barcode list from Khai, 
# format the columns to request metadata from Beverly 
# usage: metadata_setup.py <sample barcode list.csv> <output_prefix>

# read in the flowcell sample barcode list
barcode_csv = pd.read_csv(sys.argv[1], dtype=str)

# re-organize column names for Beverly
metadata_cols = barcode_csv.columns.tolist()
import_cols = {}
for col in metadata_cols:
	if 'sampleid' == col.lower():
		import_cols['SampleID'] = col
		continue
	elif 'barcode' == col.lower():
		import_cols['barcode'] = col
		continue
	elif 'pcr barcode' == col.lower():
		import_cols['pcr_barcode'] = col
		continue
	elif 'library' == col.lower():
		import_cols['library'] = col
		continue
	elif 'comments' == col.lower():
		import_cols['kn_comments'] = col
		continue

# set the columns to keep and request from Beverly
cols_to_keep = ['rfid', 'library', 'sample_id', 'project', 'pool', 'runid', 'barcode', 'pcr_barcode',
                'sex', 'coatcolor', 'sires', 'dams', 'fastq', 'kn_comments']

# add columns to the dataframe
for col in cols_to_keep:
	if col not in import_cols.keys():
		barcode_csv.insert(0, col, '')
	else:
		barcode_csv.insert(0, col, barcode_csv[import_cols[col]])

# create the sample_id by joining library and RFID
barcode_csv['sample_id'] = barcode_csv['library','rfid'].agg('_'.join, axis=1)

# save to file
barcode_csv.to_csv(sys.argv[2] + '.csv', index=False, sep=',')
