import os
import sys
import argparse
import datetime
import pandas as pd

def main(args):
    """ Main function to run all the process """
    input_path = os.path.abspath(args.input_path)
    tissue = str(args.tissue)
    output_path = os.path.abspath(args.output_path)
    # datetime.datetime.now().strftime('%y-%m-%d %a %H:%M:%S')
    datetime_now = datetime.datetime.now().strftime('%y-%m-%d')
    # Check if directories exist
    if os.path.isdir(input_path) == False:
        print(f"The input path doesn't exist: {input_path}")
        exit()
    elif os.path.isdir(output_path) == False:
        print(f"The output path doesn't exist: {output_path}")
        exit()
    else:
        print(f"Input path: {input_path}")
        print(f"Output path: {output_path}")
    output_file = os.path.join(output_path, "featureCounts_" + tissue + f"_{datetime_now}.tsv")
    dir_list = os.listdir(input_path)
    counts_list = [counts for counts in dir_list if counts.endswith(".exon.ct.short.txt")]

    # Open all tsv files and concatenate
    df_list = []
    for counts in counts_list:
        sample_name = counts.split(".")[0]
        df = pd.read_csv(os.path.join(input_path, counts), sep=" ", skiprows=2, names=["", sample_name])
        df.set_index("", inplace=True)
        df_list.append(df)
    all_dfs = pd.concat(df_list, axis=1)
    # Create final TSV table
    all_dfs.to_csv(output_file, index=True, sep="\t")

if __name__ == "__main__":
    """ Arguments for the main function """
    parser = argparse.ArgumentParser(
        description="Script to concatenate all the counts files available in one directory", usage='''Usage:
        python concat_counts.py -p /path/to/featureCounts_fibroblast/''')
    parser.add_argument("-i", "--input_path", help="Path to Counts TSV files (expecting files ending with `.exon.ct.short.txt`)", type=str, required=True)
    parser.add_argument("-t", "--tissue", help="Tissue used (eg.: fibroblast, blood, fat, muscle, etc.)", type=str, default="fibroblast")
    parser.add_argument("-o", "--output_path", help="Path to create the combined count table TSV file", type=str, default=".")
    parser.set_defaults(func=main)
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()
    args.func(args)