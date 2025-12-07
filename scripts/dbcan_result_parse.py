import pandas as pd

## data load part
# read gene_count_matrix.tsv 
gene_count_matrix = pd.read_csv('ppanggolin_result/gene_count_matrix.tsv', sep='\t')
# dbCAN all result combined 
dbcan3_result = pd.read_csv('db_can_out/overview.txt', sep='\t')

#this file is already converted into GeneID\t every CAZYme family in one gene. 
protein_family_info = pd.read_csv('dbcan_family.tsv', sep='\t')
# Rename 'Gene ID' to 'Gene' in dbcan_family.tsv header 
protein_family_info.rename(columns={'Gene ID': 'Gene'}, inplace=True)


raw_data_count = pd.merge(gene_count_matrix, dbcan3_result, left_on='Gene', right_on='Gene ID', how='outer')
raw_data_count.to_csv("dbcan_raw_gene_count_data.csv", index=False)
print("Raw data file created: dbcan_raw_gene_count_data.csv")

# Combine protein_family_info and gene_count_matrix
merged_family_count = pd.merge(protein_family_info, gene_count_matrix, on='Gene', how='left')

merged_family_count.to_csv("dbcan_geneID_HMMER_count_gene_count.csv", index=False)
print("Gene ID, HMMER, and count data file created: dbcan_geneID_HMMER_count_gene_count.csv")

# Get CAZyme family counts (considering duplicates)
# save the file 
hmmer_count_gene_count_matrix = merged_family_count.groupby('HMMER').sum().drop('Gene', axis=1, errors='ignore')
hmmer_count_gene_count_matrix.to_csv("dbcan_HMMER_count_gene_count_matrix.csv")
print("CAZyme family count file created: dbcan_HMMER_count_gene_count_matrix.csv")

# gene pa part
# gene pa part
gene_pa_matrix = pd.read_csv('ppanggolin_result/gene_presence_absence.Rtab', sep='\t')


# Create raw data by merging gene_pa_matrix and dbcan3_result
#dbcan3_result called upper part
raw_data_pa = pd.merge(gene_pa_matrix, dbcan3_result, left_on='Gene', right_on='Gene ID', how='outer')
raw_data_pa.to_csv("dbcan_raw_gene_PA_data.csv", index=False)
print("Raw data file created: dbcan_raw_gene_PA_data.csv")

# Combine protein_family_info and gene_pa_matrix
merged_family_pa = pd.merge(protein_family_info, gene_pa_matrix, on='Gene', how='left')

# Save merged data with Gene ID, HMMER, and gene PA
merged_family_pa.to_csv("dbcan_geneID_HMMER_count_gene_PA.csv", index=False)
print("Gene ID, HMMER, and PA data file created: dbcan_geneID_HMMER_count_gene_PA.csv")

# Get CAZyme family PA counts (considering duplicates)
hmmer_pa_matrix = merged_family_pa.groupby('HMMER').sum().drop('Gene', axis=1, errors='ignore')

# Save HMMER PA matrix
hmmer_pa_matrix.to_csv("dbcan_HMMER_count_gene_PA_matrix.csv")
print("CAZyme family PA file created: dbcan_HMMER_count_gene_PA_matrix.csv")

print("Processing completed.")