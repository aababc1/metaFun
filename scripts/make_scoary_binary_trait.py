import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Process TSV file and create gene traits.')
parser.add_argument('-i', '--input', required=True, help='input metadata file')
parser.add_argument('-s', '--samplecol', type=int, default=1, help='column number for sample ID (1-based)')
parser.add_argument('-m', '--metacol', type=int, required=True, help='column number for selected metadata (1-based)')

args = parser.parse_args()

#df = pd.read_csv(args.input, sep='\t')
df = pd.read_csv(args.input)

sample_col_idx = args.samplecol - 1
try:
    sample_col = df.columns[sample_col_idx]
except IndexError:
    print(f"Error: Invalid sample ID column number. Please specify a valid sample ID column number.")
    exit(1) 

# metadata column selection
meta_col_idx = args.metacol - 1
try:
    meta_col = df.columns[meta_col_idx]
except IndexError:
    print(f"Error: Invalid metadata column number. Please specify a valid column number.")
    exit(1)
# extract needed columns 
selected_df = df[[sample_col, meta_col]]

category_df = pd.get_dummies(selected_df.iloc[:, 1])  # second one is metadata column ID 
result_df = pd.concat([selected_df.iloc[:, 0], category_df], axis=1)
result_df = result_df * 1

# generate input file for scoary2  
result_df.to_csv('genome_traits.tsv', sep='\t', index=False)

