#!/usr/bin/env python

import os
import sys
import argparse
import datetime
import pandas as pd

def main(args):
    """ Main function to run all the process """
    master_file = os.path.abspath(args.master_file)
    new_sample = os.path.abspath(args.new_sample)
    output_path = os.path.abspath(args.output_path)
    # datetime.datetime.now().strftime('%y-%m-%d %a %H:%M:%S')
    datetime_now = datetime.datetime.now().strftime('%y-%m-%d')
    output_name = "_".join(os.path.basename(master_file).split("_")[:2])
    print(output_name)
    output_file = os.path.join(output_path, output_name  + f"_{datetime_now}_updated.tsv")


    # Open files to combine
    df_master = pd.read_csv(master_file, sep="\t")
    df_master.set_index("Unnamed: 0", inplace=True)
    print(df_master.head())
    df_master.index.names = ['']
    print(df_master.head())
    sample_name = os.path.basename(new_sample).split(".")[0]
    df_sample = pd.read_csv(new_sample, sep=" ", skiprows=2, names=["", sample_name])
    df_sample.set_index("", inplace=True)
    print(df_sample.head())
    
    # Merge files and save it
    all_dfs = pd.concat([df_master, df_sample], axis=1)
    all_dfs.to_csv(output_file, index=True, sep="\t")

if __name__ == "__main__":
    """ Arguments for the main function """
    parser = argparse.ArgumentParser(
        description="Script to concatenate all the counts files available in one directory", usage='''Usage:
        python add1_count.py -m /path/to/featureCounts_fibroblast.csv -n path/to/UDN208852-M-fibroblast.gene_id.exon.ct.short.txt''')
    parser.add_argument("-m", "--master_file", help="Master file with all controls samples", type=str, required=True)
    parser.add_argument("-n", "--new_sample", help="New sample count file (expecting files ending with `.exon.ct.short.txt`)", type=str, required=True)
    parser.add_argument("-o", "--output_path", help="Path to create the combined count table TSV file", type=str, default=".")
    parser.set_defaults(func=main)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)