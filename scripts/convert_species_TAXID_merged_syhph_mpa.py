import pandas as pd
import argparse
import sys
import numpy as np

def load_taxonomy_mapping(metadata_file):
    """
    
    Parameters:
        metadata_file (str): metadata path
    
    Returns:
        dict: { 's__Species Name': taxonomy_id, ... }
    """
    try:
        # sylph accession tax taxID file .
        # read metadata file 
        df = pd.read_csv(metadata_file, sep='\t', dtype=str)
    except Exception as e:
        print(f"ERROR: Error during reading genome taxID file : {e}")
        sys.exit(1)
    
    # check essential rows
    required_columns = {'accession', 'gtdb_taxonomy', 'taxonomy_id'}
    if not required_columns.issubset(df.columns):
        print(f"ERROR: File for sylph genome db taxID mapping is wrong.  Required rows: {required_columns}")
        sys.exit(1)
    
    species_taxa_to_id = {}
    
    for _, row in df.iterrows():
        gtdb_taxonomy = row['gtdb_taxonomy']
        taxonomy_id = row['taxonomy_id']
        # extract species 
        taxa_levels = gtdb_taxonomy.split(';')
        if len(taxa_levels) < 7:
            print(f"WARNING: gtdb_taxonomy error: {gtdb_taxonomy}")
            continue
        
        species_taxa = taxa_levels[6].strip()  #7th element (0-based index: 6)
        if not species_taxa.startswith('s__'):
            continue
        
        #map to dictionary
        species_taxa_to_id[species_taxa] = taxonomy_id
    
    print(f"INFO: metadata load completed. Number of all species : {len(species_taxa_to_id)}")
    return species_taxa_to_id

def process_sylph_mpa(mpa_file, taxonomy_mapping, output_file):
    """
    read sylphmpa file and keep species rows convert taxonmy to taxonomy_id
    Parameters:
        mpa_file (str): Sylph MPA input file path .
        taxonomy_mapping (dict): species_taxa_to_id mapping dictionary.
        output_file (str): final output path.
    """
    try:
        with open(mpa_file, 'r') as infile, open(output_file, 'w') as outfile:
            # read first line
            header = infile.readline().strip().split('\t')
            if len(header) < 2:
                print("ERROR: sylphmpa file is wrong.")
                sys.exit(1)
            # set taxonomy_id as first row in output file 
            new_header = ['taxonomy_id'] + header[1:]
            outfile.write('\t'.join(new_header) + '\n')
            
            line_count = 0
            mapped_count = 0
            skipped_count = 0
            
            for line in infile:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    print(f"WARNING: lines do not contain enough columns. skip line: {line.strip()}")
                    skipped_count += 1
                    continue
                
                taxonomy = parts[0].strip()
               # species level , classification 
               
                if '|s__' not in taxonomy and not taxonomy.startswith('s__'):
                    skipped_count += 1
                    continue
                
                # extract starts with 's__'
                if taxonomy.startswith('s__'):
                    species_taxa = taxonomy
                else:
                    #  '|s__'
                    species_taxa = taxonomy.split('|')[-1].strip()
                    
                # taxonomy_id mapping
                taxonomy_id = taxonomy_mapping.get(species_taxa, None)
                if taxonomy_id:
                    try:
                        values = np.array(parts[1:], dtype=np.float64)
                        values = np.where(values != 0, values / 100, 0)
                        values = np.round(values, 6)

                        str_values = [f"{v:.6f}" for v in values]
                        output_parts = [taxonomy_id] + str_values
                        outfile.write('\t'.join(output_parts) + '\n')
                        mapped_count += 1
                    except ValueError as e:
                        print(f"WARNING: Error processing numerical values: {e}")
                        skipped_count += 1
                else:
                    print(f"WARNING: taxonomy '{species_taxa}' is not at the metadata mapping file. skip line")
                    skipped_count += 1

                # # taxonomy_id mapping
                # taxonomy_id = taxonomy_mapping.get(species_taxa, None)
                # if taxonomy_id:
                #     # print lines 
                #     output_parts = [taxonomy_id] + parts[1:]
                #     outfile.write('\t'.join(output_parts) + '\n')
                #     mapped_count += 1
                # else:
                #     print(f"WARNING: taxonomy '{species_taxa}' is not at the metadata mapping file. skip line")
                #     skipped_count += 1
                
                line_count += 1
            
            print(f"INFO: total lines: {line_count}")
            print(f"INFO: mapped lines: {mapped_count}")
            print(f"INFO: skipped lines: {skipped_count}")
            print(f"INFO: output file : {output_file}")
    
    except FileNotFoundError as fe:
        print(f"ERROR: Cannot find file: {fe}")
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Error during handling file: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Sylph MPA file to  taxonomy_id.")
    parser.add_argument("-m", "--metadata", help="metadata file path (BrackentaxID_GTDBtaxa_completed.csv.cs).", type=str, required=True)
    parser.add_argument("-i", "--input-mpa", help="Sylph MPA input file path.", type=str, required=True)
    parser.add_argument("-o", "--output", help="output file path.", type=str, required=True)
    
    args = parser.parse_args()
    
    # load metadata mapping
    taxonomy_mapping = load_taxonomy_mapping(args.metadata)
    
    # Sylph MPA file handeling 
    process_sylph_mpa(args.input_mpa, taxonomy_mapping, args.output)

if __name__ == "__main__":
    main()
