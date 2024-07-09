import pandas as pd
import sys

quality_report_path = sys.argv[1]
gunc_report_path = sys.argv[2]
output_file_path = sys.argv[3]

# Load the data from the files
#quality_report = pd.read_csv('quality_report.tsv', sep='\t')
quality_report = pd.read_csv(quality_report_path, sep='\t')

#gunc_report = pd.read_csv('GUNC.progenomes_2.1.maxCSS_level.tsv', sep='\t')
gunc_report = pd.read_csv(gunc_report_path, sep='\t')
# Rename the columns for a clear key to merge
quality_report.rename(columns={'Name': 'Genome'}, inplace=True)

# Merge the dataframes on the Genome column
merged_data = pd.merge(quality_report[['Genome', 'Completeness', 'Contamination']],
                       gunc_report[['genome', 'clade_separation_score', 'contamination_portion', 
                                    'reference_representation_score', 'pass.GUNC']],
                       left_on='Genome', right_on='genome')

# Drop the duplicate 'genome' column
merged_data.drop(columns=['genome'], inplace=True)

# QS50 filter : Completeness - 5*Contamination >= 50 
# Calculate QS50  
merged_data['QS50'] = merged_data['Completeness'] - 5 * merged_data['Contamination']
merged_data['QS50.pass'] = (merged_data['QS50']) >= 50 & (merged_data['Contamination'] < 5) & (merged_data['Completeness'] >= 50) 
# Define medium and near_complete quality (following MIMAG quality, tRNA and rRNA is omitted) with GUNC pass
merged_data['medium_quality_gunc.pass'] = (merged_data['Completeness'] >= 50) & (merged_data['Contamination'] < 10)  & (merged_data['pass.GUNC'] == True)
merged_data['near_complete_gunc.pass'] = (merged_data['Completeness'] > 90) & (merged_data['Contamination'] < 5)  & (merged_data['pass.GUNC'] == True)
merged_data['medium_quality.pass'] = (merged_data['Completeness'] >= 50) & (merged_data['Contamination'] < 10) 
merged_data['near_complete.pass'] = (merged_data['Completeness'] > 90) & (merged_data['Contamination'] < 5) 


# quality pass based on checkm2 completeness contamination ,  QS50, gunc pass 
QS50_GUNC = (merged_data['Contamination'] < 5) & (merged_data['Completeness'] >= 50) & \
          (merged_data['QS50.pass']) & (merged_data['pass.GUNC'] == True)
merged_data['QS50_gunc.pass'] = QS50_GUNC

column_order = ['Genome', 'Completeness', 'Contamination','medium_quality.pass','near_complete.pass', 'medium_quality_gunc.pass', 'near_complete_gunc.pass', 'QS50', 'QS50.pass', 'pass.GUNC', 'QS50_gunc.pass']
merged_data = merged_data[column_order]

# Save the merged data to a new file
#merged_data.to_csv('combined_report.tsv', sep='\t', index=False)

merged_data.to_csv(output_file_path, sep='\t', index=False)
