import pandas as pd

# read gene_count_matrix.tsv 
gene_count_matrix = pd.read_csv('ppanggolin_result/gene_count_matrix.tsv', sep='\t')

dbcan3_result = pd.read_csv('db_can_out/overview.txt', sep='\t')

# get geneID list using dbcan3_result in gene_count matrix from ppangglin
merged_data = pd.merge(gene_count_matrix, dbcan3_result, left_on='Gene', right_on='Gene ID', how='right')

# dbCAN results are combined with individual gene table matrix
merged_data.to_csv("dbcan_gene_count_matrix.csv", index=False)

# dbcan_family information table
protein_family_info = pd.read_csv('dbcan_family.tsv', sep='\t')
protein_family_info.rename(columns={'Gene ID': 'Gene'}, inplace=True)

# combine dbcan_family and  gene_count_matrix 
merged_data = pd.merge(protein_family_info, gene_count_matrix, on='Gene')

merged_dbcan_family = merged_data.groupby(['HMMER', 'Gene']).sum()
#This will not be stored
#merged_dbcan_family.to_csv("dbcan_gene_dbcan_family.csv", index=True)


#Get CAZyme family counts
hmmer_count_matrix = merged_data.groupby(['HMMER', 'Gene']).sum().groupby('HMMER').sum()

#print result 
hmmer_count_matrix.to_csv("dbcan_hmmerfamily_count_matrix.csv")
