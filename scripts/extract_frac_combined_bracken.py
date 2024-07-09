import argparse

def extract_and_write_desired_columns(input_file_path, output_file_path):
    with open(input_file_path, 'r') as infile, open(output_file_path, 'w') as outfile:
        header_line = infile.readline().strip()
        columns = header_line.split('\t')
        frac_columns_indices = [i for i, column in enumerate(columns[3:], start=3) if column.endswith('_frac')]
        desired_column_indices = [1] + frac_columns_indices
        desired_columns = [columns[i] for i in desired_column_indices]
        outfile.write('\t'.join(desired_columns) + '\n')
        for line in infile:
            data = line.strip().split('\t')
            desired_data = [data[i] for i in desired_column_indices]
            outfile.write('\t'.join(desired_data) + '\n')

if __name__ == "__main__":
    # Set up argparse to handle command line arguments
    parser = argparse.ArgumentParser(description='Extract specific columns from a TSV file.')
    parser.add_argument('input_file', help='Path to the input TSV file')
    parser.add_argument('output_file', help='Path to the output TSV file')

    # Parse the arguments
    args = parser.parse_args()

    # Call the function with the provided arguments
    extract_and_write_desired_columns(args.input_file, args.output_file)


