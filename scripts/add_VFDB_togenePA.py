import argparse
import pandas as pd
import re 

# Set up the argument parser
parser = argparse.ArgumentParser(description='Create a heatmap with VFDB information.')
parser.add_argument('-d', '--diamond_output', required=True, help='DIAMOND output file path')
parser.add_argument('-g', '--gene_pa', required=True, help='Gene presence absence file path')
parser.add_argument('-va', '--vfdb_annot_ref', required=True, help='VFDB setB annotation reference file path')
parser.add_argument('-o', '--output', required=True, help='Output enriched gene PA file path')

args = parser.parse_args()

diamond_output = pd.read_csv(args.diamond_output, sep='\t', header=None, \
                            names=['gene_name', 'vfdb_id','pident','length','mismatch','gapopen','qstart','qend','sstart','send','evalue','bitscore'])

# diamond_output['VFG_ID'] = diamond_output['vfdb_id'].str.extract(r'(VFG\d{6})')
vfdb_annot_ref = pd.read_csv(args.vfdb_annot_ref, sep='\t', header=None, names=['VFG_ID', 'gene_ab', 'gene_description', 'VF_ID', 'VFC_ID', 'organism'])
annotated_diamond_output = pd.merge(diamond_output, vfdb_annot_ref, left_on='vfdb_id',right_on='VFG_ID', how='left')
annotated_diamond_output.to_csv("geneID_CARD_annotation.csv",index =False)


# read genePA Rtab tsv from panaroo, retain annotated genes in PA file. 
gene_pa_raw = pd.read_csv(args.gene_pa, sep='\t', header = 0) #, dtype=str)
gene_names_invfdb = diamond_output['gene_name'].unique()
#gene_names = raw_data.iloc[:, 0].unique()

#filtered_gene_pa = gene_pa_raw.index.isin(gene_names_invfdb)
#filtered_gene_pa = gene_pa_raw[gene_pa_raw['Gene'].isin(gene_names_invfdb)]

filtered_gene_pa = gene_pa_raw[gene_pa_raw['Gene'].isin(diamond_output['gene_name'])]


# Merge the VF information with gene presence/absence data based on the gene_id



merged_data = pd.merge(filtered_gene_pa, annotated_diamond_output[[ 'gene_name','VFG_ID', 'gene_description', 'VF_ID', 'VFC_ID', 'organism']], 
                       left_on='Gene', right_on='gene_name', how='left')

# Save the  gene presence and absence data to an output file
merged_data.to_csv(args.output, index=False)
# # Selecting the relevant columns for the final output
# final_columns = gene_pa.columns.tolist() + ['VF_Number']
# final_data = merged_data[final_columns]

# Save the  gene presence and absence data to an output file
