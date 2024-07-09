import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Process panaroo gene presence_absence.csv ")
parser.add_argument("-i", "--input", required=True, help="Input file path for gene presence-absence data.")
parser.add_argument("-o", "--output", required=True, help="Output file path for processed data.")

args = parser.parse_args()



# Load the gene presence-absence file with only necessary columns
gene_pa_df = pd.read_csv(args.input, usecols=lambda column: column not in ['Non-unique Gene name', 'Annotation'])

# Load the KO list file and drop any missing KO values
ko_list_df = pd.read_csv('kofamscan_result', sep='\t', header=None, names=['Gene', 'KO']).dropna()

# Convert Gene columns to string to ensure proper matching
ko_list_df['Gene'] = ko_list_df['Gene'].astype(str)

# Create a dictionary for fast KO lookup
ko_dict = ko_list_df.set_index('Gene')['KO'].to_dict()

# Extract the non-empty gene entries from gene_presence_absence.csv
non_empty_genes = gene_pa_df[gene_pa_df['Gene'].notna()]

# Map the genes to their KOs using the dictionary
non_empty_genes['KO'] = non_empty_genes['Gene'].map(ko_dict)

# Now, we need to pivot the table to have one row per gene per genome
# First, let's melt the non_empty_genes dataframe
melted = non_empty_genes.melt(id_vars=['Gene', 'KO'], var_name='Genome', value_name='GeneValue')

# We only want to keep rows where GeneValue is not NaN (meaning the gene is present in the genome)
filtered = melted[melted['GeneValue'].notna()]

# Save the final DataFrame to a CSV file
filtered[['Genome', 'GeneValue', 'KO']].to_csv(args.output, index=False)

