import pandas as pd
import sys

# read in the metadata
original_metadata = pd.read_csv(sys.argv[1], dtype=str)

# filter for HS rats only
original_metadata['strain'] = original_metadata['strain'].replace('Heterogenous stock', 'Heterogeneous stock')
original_metadata = original_metadata[original_metadata['strain'] == 'Heterogeneous stock'].reset_index(drop=True)

# set data types
original_metadata['pcr_barcode'] = pd.to_numeric(original_metadata['pcr_barcode'], errors='coerce').astype(int)

# re-organize column names for Fgbio demultiplex
metadata_cols = original_metadata.columns.tolist()
import_cols = {}
for col in metadata_cols:
	if 'rfid' == col.lower():
		import_cols['Sample_Name'] = col
		continue
	elif 'library_name' == col.lower():
		import_cols['Library_ID'] = col
		continue
	elif 'project_name' == col.lower():
		import_cols['Sample_Project'] = col
		continue
	elif 'barcode' == col.lower():
		import_cols['Sample_Barcode'] = col
		continue
	# elif 'flowcell_id' == col.lower():
	# 	import_cols['full_run_id'] = col
	# 	continue
	elif 'fastq_files' == col.lower():
		import_cols['Fastq_Files'] = col
		continue
	elif 'pcr_barcode' == col.lower():
		import_cols['PCR_Barcode'] = col
		continue
for new_col in ['Fastq_Files', 'Sample_Barcode', 'Sample_Project', 'Library_ID', 'Sample_Name']:
	if new_col not in import_cols.keys():
		original_metadata.insert(0, new_col, '')
		print("ERROR: " + new_col + " doesn't exist")
	else:
		original_metadata.insert(0, new_col, original_metadata[import_cols[new_col]])

temp_sample_ID = original_metadata[[import_cols["Library_ID"], import_cols["Sample_Name"]]].agg('_'.join, axis=1)
original_metadata.insert(0, "Sample_ID", temp_sample_ID)

# save it to sample_sheet.csv
original_metadata.to_csv(sys.argv[2]+"/sample_sheet.csv", index=False, sep=',')

# seperate samples based on library to different SampleSheet.csv for demultiplex
for unique in original_metadata['Library_ID'].unique():
	temp_metadata = original_metadata[original_metadata['Library_ID'] == unique].reset_index(drop=True)
	if 'PCR_Barcode' not in import_cols.keys():
		pcr_barcode = "NONE"
		print("ERROR: pcr_barcode doesn't exist")
	else:
		pcr_barcode = str(temp_metadata.iloc[0][import_cols['PCR_Barcode']])
		# pcr_barcode = temp_metadata['PCR_Barcode'][0].astype(int)
	# if 'full_run_id' not in import_cols.keys():
	# 	full_run_id = "NONE"
	# 	print("ERROR: runid doesn't exist")
	# else:
	# 	# full_run_id = str(temp_metadata.iloc[0][import_cols['full_run_id']])
	# 	full_run_id = str(temp_metadata[['full_run_id']][0])
	# temp_metadata.to_csv(sys.argv[2]+"/SampleSheet_"+pcr_barcode+"_"+unique+"_"+full_run_id+".csv", index=False, sep=',')
		flowcell_id = str(temp_metadata['flowcell_id'][0])
	temp_metadata.to_csv(sys.argv[2]+"/SampleSheet_"+pcr_barcode+"_"+unique+"_"+flowcell_id+".csv", index=False, sep=',')
