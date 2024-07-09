import pandas as pd
import argparse
import re
import sys

def main(genome_metadata_path, metagenome_metadata_path, output_path,accession_column):
    # load metadata
    metadata_metagenome = pd.read_csv(metagenome_metadata_path)
    metadata_genome = pd.read_csv(genome_metadata_path)

    # Generate  file name prefix accession (that should be in metagenome metadata)
    prefixes = ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global', 'self']
    #prefix_pattern = '|'.join(prefixes)
    prefix_pattern = '|'.join([f"_{p}" for p in prefixes])

    #metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(r'(_SB2.*|_MB2.*)', '', x))
    #first_column_metagenome = metadata_metagenome.columns[0]
    #metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(rf"(?:{prefix_pattern})_SB2.*", '', x))
    metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(rf"(?:(?:{prefix_pattern})_SB2|_MB2).*", '', x))

    genome_accessions = set(metadata_genome['Analysis_accession'])
    metagenome_accessions = set(metadata_metagenome.iloc[:, int(accession_column) - 1])
    if not genome_accessions.issubset(metagenome_accessions):
        print("Error: The unique values in the 'Analysis_accession' column of the genome metadata do not match the unique values in the specified accession column of the metagenome metadata.")
        print("Please specify the correct accession column using the -a or --accession_column argument.")
        sys.exit(1)


    # Merge the two dataframes
    #merged_metadata = pd.merge(metadata_genome,metadata_metagenome,  left_on='Analysis_accession' , right_on=first_column_metagenome)
    #merged_metadata = pd.merge(metadata_genome, metadata_metagenome, left_on='Analysis_accession', right_on=accession_column)
    merged_metadata = pd.merge(metadata_genome, metadata_metagenome, left_on='Analysis_accession', right_on=metadata_metagenome.columns[int(accession_column) - 1],how = 'left')
    
    
    # Check column is right 
    if merged_metadata['Analysis_accession'].isnull().any():
        print("Warning: The specified accession column does not match the 'Analysis_accession' column in the genome metadata.")
        print("Please specify the correct accession column using the -a or --accession_column argument.")

    # Save the merged metadata
    merged_metadata.to_csv(output_path, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge two metadata files.")
    parser.add_argument("-g", "--genome_metadata", required=True, help="Path to the genome metadata file")
    parser.add_argument("-m", "--metagenome_metadata", required=True, help="Path to the metagenome metadata file")
    parser.add_argument("-o", "--combined_metadata", required=True, help="Output file path for the combined metadata")
#    parser.add_argument("-a", "--accession_column", required=True, help="Column name for accession in metagenome metadata")
    parser.add_argument("-a", "--accession_column", default='1', help="Column number for accession in metagenome metadata (1-based index, default: 1)")
    

    args = parser.parse_args()
    #main(args.genome_metadata, args.metagenome_metadata, args.combined_metadata)
    main(args.genome_metadata, args.metagenome_metadata, args.combined_metadata, args.accession_column)
