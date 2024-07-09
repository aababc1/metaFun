import pandas as pd
import argparse
import re

def main(genome_metadata_path, metagenome_metadata_path, output_path):
    # load metadata
    metadata_metagenome = pd.read_csv(metagenome_metadata_path)
    metadata_genome = pd.read_csv(genome_metadata_path)

    # Generate  file name prefix accession (that should be in metagenome metadata)
    metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(r'(metabat2_|semibin2_)', '', x))
    metadata_genome['Analysis_accession'] = metadata_genome['Analysis_accession'].apply(lambda x: re.sub(r'(_SB2.*|_MB2.*)', '', x))
    first_column_metagenome = metadata_metagenome.columns[0]

    # Merge the two dataframes
    merged_metadata = pd.merge(metadata_genome,metadata_metagenome,  left_on='Analysis_accession' , right_on=first_column_metagenome)

    # Save the merged metadata
    merged_metadata.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge two metadata files.")
    parser.add_argument("-g", "--genome_metadata", required=True, help="Path to the genome metadata file")
    parser.add_argument("-m", "--metagenome_metadata", required=True, help="Path to the metagenome metadata file")
    parser.add_argument("-o", "--combined_metadata", required=True, help="Output file path for the combined metadata")

    args = parser.parse_args()
    main(args.genome_metadata, args.metagenome_metadata, args.combined_metadata)
