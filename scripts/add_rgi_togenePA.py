import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Set up the argument parser
parser = argparse.ArgumentParser(description='Create a heatmap with additional ARO information.')
parser.add_argument('-i', '--input', required=True, help='Input raw data file path')
parser.add_argument('-o', '--output', required=True, help='Output heatmap file path')
parser.add_argument('-r', '--aro_ref', required=True, help='ARO reference file path')
parser.add_argument('-gpa', '--gene_pa', required=True, help='Gene presence absence file path')

# Parse the arguments
args = parser.parse_args()

# Load the raw data, gene presence/absence data, and ARO reference data
raw_data = pd.read_csv(args.input, sep='\t')
gene_pa = pd.read_csv(args.gene_pa, index_col=0,sep='\t') #, set='\t'
aro_ref = pd.read_csv(args.aro_ref, sep='\t')

# only select one of duplicated ARO_accession information 
aro_ref = aro_ref.drop_duplicates(subset=['ARO Accession'])



# Filter gene_pa using the first column (gene names) of raw_data
gene_names = raw_data['ORF_ID'].unique()
#gene_names = raw_data.iloc[:, 0].unique()
filtered_gene_pa = gene_pa.loc[gene_pa.index.isin(gene_names)]


# Map ARO to 'Resistance Mechanism' and 'CARD Short Name'
raw_data['Resistance Mechanism'] = raw_data['ARO'].map(aro_ref.set_index('ARO Accession')['Resistance Mechanism'])
raw_data['CARD Short Name'] = raw_data['ARO'].map(aro_ref.set_index('ARO Accession')['CARD Short Name'])
raw_data['Drug Class'] = raw_data['ARO'].map(aro_ref.set_index('ARO Accession')['Drug Class'])  

# Merge this information with filtered_gene_pa
#selected_columns = ['ORF_ID', 'Resistance Mechanism', 'CARD Short Name']

raw_data.set_index('ORF_ID', inplace=True)


final_data = filtered_gene_pa.merge(raw_data[['Resistance Mechanism', 'CARD Short Name','Drug Class']], 
                                    left_index=True, right_index=True, how='left')


final_data.to_csv(args.output)
# Pivot final_data for heatmap
pivot_table = final_data.pivot_table(index=final_data.index, values=filtered_gene_pa.columns, aggfunc='first')

#pivot_table = final_data.pivot_table(index=raw_data.columns[0], values=filtered_gene_pa.columns, aggfunc='first')

# Plot the heatmap
#plt.figure(figsize=(12, 8))
#sns.heatmap(pivot_table, cmap='viridis', cbar_kws={'label': 'Presence/Absence'})
#plt.show()

# Optionally, write the final_data to a file
final_data.to_csv(args.output)

