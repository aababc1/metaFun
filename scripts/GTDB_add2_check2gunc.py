import pandas as pd
import sys

quality_tmp_tsv = sys.argv[1]
gtdbtk_report_path = sys.argv[2]
output_file_path = sys.argv[3]

combined = pd.read_csv(quality_tmp_tsv, sep='\t')
gtdbtk = pd.read_csv(gtdbtk_report_path, sep='\t')

gtdbtk = gtdbtk.rename(columns={'user_genome':'Genome'})
merged = pd.merge(combined, gtdbtk, on='Genome', how='left')

merged.to_csv(output_file_path, index=False)