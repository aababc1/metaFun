import pandas as pd
import argparse
import sys

def main(genome_metadata_path, user_metadata_path, output_path, first_column_metagenome):
    # Load user metadata
    if user_metadata_path.endswith('.tsv'):
        metadata_user = pd.read_csv(user_metadata_path, sep='\t')
        user_metadata_path = user_metadata_path.replace('.tsv', '.csv')
        metadata_user.to_csv(user_metadata_path, index=False)
    else:
        metadata_user = pd.read_csv(user_metadata_path)

    # Load genome metadata
    metadata_genome = pd.read_csv(genome_metadata_path)
    
    # Determine if first_column_metagenome is a number (1-based index) or a column name
    if first_column_metagenome.isdigit():
        column_index = int(first_column_metagenome) - 1  # Convert to 0-based index
        if column_index >= len(metadata_user.columns) or column_index < 0:
            raise ValueError(f"Column index {first_column_metagenome} is out of bounds.")
        column_name = metadata_user.columns[column_index]
    else:
        column_name = first_column_metagenome
        if column_name not in metadata_user.columns:
            raise ValueError(f"The specified column '{column_name}' does not exist in the user metadata.")
    #treat overlapping header names 
    overlapping_columns = [col for col in metadata_user.columns if col in metadata_genome.columns and col != column_name]
    for col in overlapping_columns:
        metadata_user.rename(columns={col: f"user_{col}"}, inplace=True)

    # Merge the two dataframes
    merged_metadata = pd.merge(metadata_genome, metadata_user, left_on='Genome', right_on=column_name)

    # Check if the number of rows matches after merging
    if len(merged_metadata) != len(metadata_genome):
        raise ValueError("The number of rows in the merged metadata does not match the genome metadata.")

    # Save the merged metadata
    merged_metadata.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge two metadata files.")
    parser.add_argument("-g", "--genome_metadata", default='quality_taxonomy_combined.csv', help="Path to the genome metadata file (default: 'quality_taxonomy_combined.csv')")
    parser.add_argument("-u", "--user_metadata", required=True, help="Path to the metagenome metadata file")
    parser.add_argument("-o", "--combined_metadata", default='combined_metadata.csv', help="Output file path for the combined metadata (default: 'combined_metadata.csv')")
    parser.add_argument("-c", "--column", required=True, help="Column index or name in the user metadata to merge on")
    args = parser.parse_args()
    
    main(args.genome_metadata, args.user_metadata, args.combined_metadata, args.column)
