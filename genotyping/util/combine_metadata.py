#!/usr/win/env python3

# This script combines metadata from previous and current flowcells for a genotyping round

# Usage: python3 combine_metadata.py combine_metadata.py -l genotyping_log -p old_metadata.csv -m metadata_list 
#                                                        -g ref_genome -o output/dir/name

import os
import sys
import re
import numpy as np
import pandas as pd
from optparse import OptionParser

def help():
    print("====== make_genotype_log.py =====")
    print("Produce a genotyping log with summary stats and QC results")
    print("-l <genotyping log>       the genotyping log for the most recent genotyping round")
    print("-p <previous metadata>    the csv file with all metadata for the most recent genotyping round")
    print("-m <new metadata>         text file listing paths for sample sheets for each flowcell new to the current round")
    print("-g <reference genome>     the base name of the reference genome")
    print("-o <output prefix>        the output prefix for all output files")
    print("Usage: combine_metadata.py -l genotyping_log -p old_metadata.csv -m metadata_list -g ref_genome -o output/dir/name")
    sys.exit()

if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-l', type="string", nargs=1, dest="genotyping_log", help="<genotyping log from the previous round>")
    parser.add_option('-p', type="string", nargs=1, dest="prev_metadata", help="<metadata csv from the previous round>")
    parser.add_option('-m', type="string", nargs=1, dest="new_metadata", help="<list of paths to new sample sheets>")
    parser.add_option('-g', type="string", nargs=1, dest="reference_genome", help="<base name of reference genome>")
    parser.add_option('-o', type="string", nargs=1, dest="output_prefix", help="<output prefix>")
    options, args = parser.parse_args()
        
    
    if options.genotyping_log != None:
        geno_log = options.genotyping_log
    else:
        raise "Please provide a genotyping log"
    if options.new_metadata != None:
        new_md = options.new_metadata
    else:
        raise "Please provide a list of sample sheets for the current genotyping round"
    if options.prev_metadata != None:
        prev_md = options.prev_metadata
    else:
        raise "Please provide the metadata file from the previous genotyping round"
 
    if options.reference_genome != None:
        ref_gen = options.reference_genome
    else:
        raise "Please provide the path to the reference genome"
    if options.output_prefix != None:
        output_prefix = options.output_prefix
    else:
        raise "Please provide the current genotyping round"


    geno_log = pd.read_csv(geno_log)
    prev_md = pd.read_csv(prev_md)
 
    # get passing samples from previous round
    good_sequences = ['analysis', 'keep']
    use_samples = geno_log[geno_log['sample_use'].isin(good_sequences)]['sample_id']

    # add sample id to previous metadata
    prev_md['sample_id'] = prev_md['library_name']+ '_' + prev_md['rfid']
    md_cols = ['sample_id', 'rfid', 'library_name', 'project_name', 'barcode', 'pcr_barcode', 'runid',  'sex', 
            'coatcolor', 'dams', 'sires', 'organism', 'strain', 'fastq_files', 'bam']
    prev_md = prev_md[md_cols]

    # keep only previous samples that passed QC
    prev_use = prev_md[prev_md['sample_id'].isin(use_samples)]

    # function to find the bam file associated with a sample ID
    def find_bam_file(sample_name):
        for bam in bamfiles:
                if sample_name in bam:
                    bam_out = f'{bam_dir}/{bam}'
                    return bam_out
                    break

    # list of metadata columns to keep/concatenate
    md_columns = ['Sample_ID', 'rfid', 'library_name', 'project_name', 'barcode', 'pcr_barcode',
                'runid', 'sex', 'coatcolor', 'dams', 'sires', 'organism', 'strain', 'fastq_files', 'bam']

    with open(new_md, 'r') as file:
        md_paths = file.readlines()

    # emtpy list for old metadata
    new_md_list = []
    new_bams = []

    for path in md_paths:

        md_path = path.strip() 
        df = pd.read_csv(md_path)
        flowcell_dir = md_path.split('/demux')[0]
        bam_dir = f'{flowcell_dir}/{ref_gen}/bams'
        
        # check and insert missing columns with NaN values
        for col in md_columns:
            if col not in df.columns:
                df[col] = pd.Series([None] * len(df), name=col)

        # reorder columns to allow concatenation
        df = df[md_columns]
        df.rename({'Sample_ID':'sample_id'},axis=1,inplace=True)

        # list all bam files in the flowcell directory
        bamfiles = [file for file in os.listdir(bam_dir) if file.endswith('.bam')]

        # add bam files to metadata
        df['bam'] = df['sample_id'].apply(find_bam_file)

        # save new bam files
        for sample in df['sample_id'].tolist(): 
            new_bams.append(find_bam_file(sample))

        new_md_list.append(df)

    # concatenate all new metradata into a single dataframe
    all_new_md = pd.concat(new_md_list, ignore_index=True)
    new_bams = pd.DataFrame({'bam' : new_bams})

    # concatenate all metadata
    all_md = pd.concat([prev_use, all_new_md], ignore_index=True)
    all_md['pcr_barcode'] = all_md['pcr_barcode'].fillna(0).astype(int)
    df['pcr_barcode'] = df['pcr_barcode'].replace(0, np.nan)


    # save combined metadata to file
    print(f'Writing metadata to {output_prefix}_metadata.csv')
    all_md.to_csv(f'{output_prefix}_metadata.csv', index=False)

    # save bam lists to file
    print(f'Writing bamlist to {output_prefix}_bamlist')
    all_md['bam'].to_csv(f'{output_prefix}_bamlist', index=False, header=False)
    print(f'Writing new bams to {output_prefix}_newbams')
    new_bams['bam'].to_csv(f'{output_prefix}_newbams', index=False, header=False)

    # save sample IDs to file
    print(f'Writing sample list to {output_prefix}_sample_ids')
    all_md['sample_id'].to_csv(f'{output_prefix}_sample_ids', index=False, header=False)


    # check that output files are formatted as desired
    print('Checking inputs...')
    idfile = f'{output_prefix}_sample_ids'
    bamlist = f'{output_prefix}_bamlist'

    with open(idfile, "r") as sample_file, open(bamlist, "r") as path_file:

        sample_ids = [line.strip() for line in sample_file]
        bam_files = [line.strip() for line in path_file]

    # ensure that all files are present in their stated directories
    missing_bams = []
    for bam_path in bam_files:
        if not os.path.exists(bam_path):
            missing_bams.append(bam_path)
    
    if missing_bams:
        print('Error: The following .bam files do were not found:')
        for path in missing_bams:
            print(path)
        print(f'Re-map the above {len(missing_bams)} samples before proceeding')
        sys.exit()
    else:
        print('All bams found')

    # ensure that all bam files are indexed
    needs_index = []
    for bam_path in bam_files:
        index_file = bam_path[:-1] + 'i'
        if not os.path.exists(index_file):
            needs_index.append(bam_path)

    if needs_index:
        print('Error: The following bam files lack a .bai index and require indexing:')
        for path in needs_index:
            print(path)
        print(f'Index the above {len(needs_index)} samples before proceeding')
        sys.exit()
    else:
        print('All bams are indexed')

    # check if the number of sample IDs and file paths match
    if len(sample_ids) != len(bam_files):
        print("Mismatch: The number of sample IDs and file paths is different")
    else:
        # check if each sample ID corresponds to the expected file path
        for i, (sample_id, file_path) in enumerate(zip(sample_ids, bam_files)):
            
            if sample_id not in file_path:
            
                print(f"Mismatch on line {i+1}: Sample ID {sample_id} does not correspond to the bam file {file_path}")
                break
        else:
            print("All samples/files correspond. You can now proceed with genotyping!")
