#!/usr/bin/env python3

import pandas as pd
import argparse
import os
import sys

def main(input_file, output_file):
    # Check if the output file already exists
    if os.path.isfile(output_file):
        print(f"Error: The output file '{output_file}' already exists. Exiting without running the script.")
        sys.exit(1)

    # Load the DataFrame from the CSV file
    df = pd.read_csv(input_file)

    # Drop rows with missing 'KO' values
    df = df.dropna(subset=['KO'])

    # Drop duplicates based on 'KO' and 'Genome' to keep only unique KO entries for each Genome
    df = df.drop_duplicates(subset=['KO', 'Genome'])

    # Pivot the DataFrame to create a matrix format with unique KO values only
    ko_matrix = df.pivot(index='KO', columns='Genome', values='KO')

    # Replace NaN with 'NA' to indicate absence
    ko_matrix.fillna('NA', inplace=True)

    # Save the resulting matrix to a CSV file
    ko_matrix.to_csv(output_file, index=True)

if __name__ == "__main__":
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Convert KO and Genome data to a matrix with unique KO values.')
    parser.add_argument('-i', '--input', required=True, help='Input CSV file path.')
    parser.add_argument('-o', '--output', required=True, help='Output CSV file path.')

    # Parse arguments
    args = parser.parse_args()

    # Run the main function with the provided arguments
    main(args.input, args.output)

