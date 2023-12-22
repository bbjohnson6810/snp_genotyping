#!/usr/win/env python3

# This script combines mapping statistics from previous and current genotyping rounds

# Usage: python3 concat_mapstats.py .py -m new_mapstats.csv -p previous_mapstats.csv -l prev_genotyping_log.csv -o output.csv

import os
import sys
import numpy as np
import pandas as pd
from optparse import OptionParser

def help():
    print("====== concat_mapstats.py =====")
    print("Produce a genotyping log with summary stats and QC results")
    print("-m <new mapstats>         the csv file with mapping statistics for all samples new to the current genotyping round")
    print("-p <previous mapstats>    the csv file with all mapping statistics for the most recent genotyping round")
    print("-l <genotyping log>       the genotyping log for the most recent genotyping round")
    print("-o <output path>        the complete output path")
    print("Usage: concat_mapstats.py --m new_mapstats.csv -p previous_mapstats.csv -l prev_genotyping_log.csv -o output.csv")
    sys.exit()

if __name__=="__main__":
    usage = "usage: %prog [options]"
    parser = OptionParser(usage)
    parser.add_option('-m', type="string", nargs=1, dest="new_mapstats", help="<new mapstats csv for the current round>")
    parser.add_option('-p', type="string", nargs=1, dest="prev_mapstats", help="<all mapstats csv from the previous round>")
    parser.add_option('-l', type="string", nargs=1, dest="genotyping_log", help="<genotyping log from the previous round>")
    parser.add_option('-o', type="string", nargs=1, dest="output_path", help="<output path>")
    options, args = parser.parse_args()
    
    if options.genotyping_log != None:
        geno_log = options.genotyping_log
    else:
        raise "Please provide a genotyping log"
    if options.new_mapstats != None:
        new_mapstats = options.new_mapstats
    else:
        raise "Please provide a csv of new mapping statistics"
    if options.prev_mapstats != None:
        prev_mapstats = options.prev_mapstats
    else:
        raise "Please provide mapping statistics from the previous genotyping round"
    if options.output_path != None:
        out_path = options.output_path
    else:
        raise "Please provide a path to the output file"

    geno_log = pd.read_csv(geno_log)
    prev_mapstats = pd.read_csv(prev_mapstats)
    new_mapstats = pd.read_csv(new_mapstats)
 
    # get passing samples from previous round
    good_sequences = ['analysis', 'keep']
    use_samples = geno_log[geno_log['sample_use'].isin(good_sequences)]['sample_id']

    # keep only previous samples that passed QC
    prev_use = prev_mapstats[prev_mapstats['sample_id'].isin(use_samples)]

    # concatenate previous & new mapping statistics
    all_mapstats_list = [prev_use, new_mapstats]
    all_mapstats = pd.concat(all_mapstats_list, ignore_index=True)

    all_mapstats = all_mapstats.drop_duplicates()

    # output to file
    print(f'Writing mapping statistics to {out_path}')
    all_mapstats.to_csv(out_path, index=False)
