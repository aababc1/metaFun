import pandas as pd
import argparse
import re
import sys

def detect_delimiter(file_path):
    with open(file_path, 'r') as f:
        first_line = f.readline()
        if '\t' in first_line:
            return '\t'
        elif ',' in first_line:
            return ','
        else:
            return ','  # default to comma if unclear

def read_metadata_file(file_path):
    """Read metadata file with automatic delimiter detection"""
    delimiter = detect_delimiter(file_path)
    try:
        return pd.read_csv(file_path, sep=delimiter)
    except Exception as e:
        print(f"Error reading file {file_path} with delimiter '{delimiter}': {str(e)}")
        sys.exit(1)

def extract_accession(genome_name):
    """
    다양한 형식의 genome 이름에서 accession을 추출하는 함수
    """
    # 기본 패턴들
    env_prefixes = ['human_gut','dog_gut','ocean','soil','cat_gut','human_oral',
                'mouse_gut','pig_gut','built_environment','wastewater',
                'chicken_caecum','global','self']

    base_patterns = [
        r'^(.+?)(?:_SB2|_MB2)',
        # NCBI (예: GCA_123456789.1)
        r'(GCA_\d+\.\d+|GCF_\d+\.\d+)',
        # number
        r'^([A-Za-z0-9]+)_'
    ]
    
    # Extract the base part first
    base_accession = None
    for pattern in base_patterns:
        match = re.search(pattern, genome_name)
        if match:
            base_accession = match.group(1)
            break
    
    if not base_accession:
        return genome_name
    
    for prefix in env_prefixes:
        if base_accession.endswith(f"_{prefix}"):
            return base_accession[:-len(f"_{prefix}")]
        elif f"_{prefix}_" in base_accession:
            return base_accession.split(f"_{prefix}_")[0]
    
    return base_accession    


    # patterns = [
    #     r'^(.+?)(?:_SB2|_MB2)',
    #     # NCBI (예: GCA_123456789.1)
    #     r'(GCA_\d+\.\d+|GCF_\d+\.\d+)',
    #     # number
    #     r'^([A-Za-z0-9]+)_'#,
    #     # 
    #     #r'^([^_]+)'
    # ]
    
    # for pattern in patterns:
    #     match = re.search(pattern, genome_name)
    #     if match:
    #         return match.group(1)
    
    # # 아무 패턴도 매치되지 않으면 원본 반환
    # return genome_name


def main(genome_metadata_path, metagenome_metadata_path, output_path,accession_column,data_type):
    # load metadata
#    metadata_metagenome = pd.read_csv(metagenome_metadata_path)
 #   metadata_genome = pd.read_csv(genome_metadata_path)
    try:
        metadata_metagenome = read_metadata_file(metagenome_metadata_path)
        metadata_genome = read_metadata_file(genome_metadata_path)
    except Exception as e:
        print(f"Error reading metadata files: {str(e)}")
        sys.exit(1)


    metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(extract_accession)
    # accession 

    print("\nAccession extraction examples:")
    sample_genomes = metadata_genome['Genome'].head()
    for genome, accession in zip(sample_genomes, metadata_genome['Analysis_accession'].head()):
        print(f"Original: {genome} -> Extracted: {accession}")
    print("=" * 50)


    # Generate  file name prefix accession (that should be in metagenome metadata)
    prefixes = ['human_gut', 'dog_gut', 'ocean', 'soil', 'cat_gut', 'human_oral', 'mouse_gut', 'pig_gut', 'built_environment', 'wastewater', 'chicken_caecum', 'global', 'self']
    #prefix_pattern = '|'.join(prefixes)
    prefix_pattern = '|'.join([f"_{p}" for p in prefixes])

    #metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(r'(_SB2.*|_MB2.*)', '', x))
    #first_column_metagenome = metadata_metagenome.columns[0]
    #metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(rf"(?:{prefix_pattern})_SB2.*", '', x))
    #metadata_genome['Analysis_accession'] = metadata_genome['Genome'].apply(lambda x: re.sub(rf"(?:(?:{prefix_pattern})_SB2|_MB2).*", '', x))

    #
    #print(f"\nFile format detection:")
    #print(f"Genome metadata delimiter: {'tab' if detect_delimiter(genome_metadata_path) == '\t' else 'comma'}")
    #print(f"Metagenome metadata delimiter: {'tab' if detect_delimiter(metagenome_metadata_path) == '\t' else 'comma'}\n")

    genome_accessions = set(metadata_genome['Analysis_accession'])
    metagenome_accessions = set(metadata_metagenome.iloc[:, int(accession_column) - 1])

    #print(f"Total genome accessions: {len(genome_accessions)}")
   # print(f"Total metagenome accessions: {len(metagenome_accessions)}")



    # Merge the two dataframes
    #merged_metadata = pd.merge(metadata_genome,metadata_metagenome,  left_on='Analysis_accession' , right_on=first_column_metagenome)
    
    merged_metadata = pd.merge(metadata_genome,
                                metadata_metagenome, 
                                left_on='Analysis_accession', 
                                right_on=metadata_metagenome.columns[int(accession_column) - 1],
                                how = 'inner')
    merged_metadata['data_type'] = data_type


    metagenome_accession_column = metadata_metagenome.columns[int(accession_column) - 1]


    total_genome = len(metadata_genome)
    total_metagenome = len(metadata_metagenome)
    used_genome = len(merged_metadata)
    used_metagenome = len(set(merged_metadata[metagenome_accession_column]))
    matched_rows = len(merged_metadata)


    print(f"Look at the file combined_metadata_quality_taxonomy.csv on your Nextflow launch directory.\n")
    print(f"Default location of output files are on results/metagenome/BIN_ASSESSMENT.\n")

    print(f"Total rows in genome metadata: {total_genome}")
    print(f"Total rows in metagenome metadata: {total_metagenome}")
    print(f"Unique genome accessions used in merged data: {used_genome}")
    print(f"Unique metagenome accessions used in merged data: {used_metagenome}")
    print(f"Total matched rows: {matched_rows}")
    print(f"Unused rows in genome metadata: {total_genome - used_genome}")
    print(f"Unused rows in metagenome metadata: {total_metagenome - used_metagenome}")

    unused_genome = set(metadata_genome['Analysis_accession']) - set(merged_metadata['Analysis_accession'])
    unused_metagenome = set(metadata_metagenome[metagenome_accession_column]) - set(merged_metadata[metagenome_accession_column])

    print("\nSample of unused accessions in genome metadata:")
    print(list(unused_genome))  
    print("\nSample of unused accessions in metagenome metadata:")
    print(list(unused_metagenome))  



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
    parser.add_argument("-d", "--data_type", default='MAG', help="Value to be added in the data_type column ")
#required=True,
    

    args = parser.parse_args()
    #main(args.genome_metadata, args.metagenome_metadata, args.combined_metadata, args.accession_column)
    main(args.genome_metadata, args.metagenome_metadata, args.combined_metadata, args.accession_column, args.data_type)

