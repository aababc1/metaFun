import os 
import argparse
# Load gene count data and metadata

parser = argparse.ArgumentParser(description='PCoA plot with PERMANOVA')
parser.add_argument('--cpus', type=int, default=8, help='Number of CPU threads to use (default: 8)')
args = parser.parse_args()

os.environ["OMP_NUM_THREADS"] = str(args.cpus)
os.environ["MKL_NUM_THREADS"] = str(args.cpus)  # Intel MKL
os.environ["NUMEXPR_NUM_THREADS"] = str(args.cpus)  # Numexpr
os.environ["NUMBA_NUM_THREADS"] = str(args.cpus)  # Numba
os.environ["NUMPY_NUM_THREADS"] = str(args.cpus)  # NumPy directly (if applicable)
import pandas as pd
import numpy as np
from skbio.diversity import beta_diversity
from skbio.stats.ordination import pcoa
from skbio.stats.distance import permanova
import plotly.graph_objects as go
from plotly.offline import plot
import plotly.express as px


gene_counts = pd.read_csv('ppanggolin_result/gene_count_matrix.tsv',sep='\t', index_col=0)
metadata = pd.read_csv('selected_metadata.csv', index_col=0)



# Calculate Bray-Curtis distance matrix and PCoA
bray_dm = beta_diversity("braycurtis", gene_counts.T)
bray_pcoa = pcoa(bray_dm)

# Prepare PCoA results for plotting
pcoa_results = pd.DataFrame(bray_pcoa.samples.iloc[:, :2], columns=['PC1', 'PC2'])
pcoa_results['sample_id'] = pcoa_results.index
merged_data = pcoa_results.merge(metadata, left_on='sample_id', right_index=True)

# Select columns for PERMANOVA
permanova_columns = metadata.columns[np.r_[0, 1, 14:len(metadata.columns)]]

# Dictionary to store PERMANOVA results and cleaned data
permanova_results = {}
cleaned_data_dict = {}

# Clean metadata and calculate PERMANOVA results
for col in permanova_columns:
    # Remove rows with NA values in the current column
    cleaned_data = merged_data.dropna(subset=[col])
    cleaned_data_dict[col] = cleaned_data

    # Report the number of NA values removed
    na_count = len(merged_data) - len(cleaned_data)
    print(f"Column '{col}' had {na_count} NA values removed.")

    # Check if there are enough samples left after NA removal
    if cleaned_data[col].nunique() > 1:
        # Ensure distance matrix and cleaned_data have matching IDs
        valid_ids = cleaned_data['sample_id'].values
        bray_dm_subset = bray_dm.filter(valid_ids)

        # Calculate PERMANOVA
        try:
            result = permanova(bray_dm_subset, cleaned_data, column=col, permutations=999)
            print(result)
            r2 = 1 - 1 / (1 + result[4] * result[3] / (result[2] - result[3] - 1))
            result['R²'] = r2
            permanova_results[col] = result
        except ValueError as e:
            print(f"Skipping {col} due to error: {e}")
    else:
        print(f"Skipping {col} due to insufficient unique values after NA removal.")

# Create dropdown menu for PERMANOVA results
def create_hover_text(df):
    hover_text = df.apply(lambda row: '<br>'.join([f"{col}: {row[col]}" for col in df.columns]), axis=1)
    return hover_text

permanova_dropdown = go.layout.Updatemenu(
    buttons=list([
        dict(
            args=[{
                'marker': {
                    'color': cleaned_data_dict[col][col].astype('category').cat.codes if cleaned_data_dict[col][col].dtype == 'category' or cleaned_data_dict[col][col].dtype == 'object' else cleaned_data_dict[col][col],
                    'colorscale': 'Viridis' if cleaned_data_dict[col][col].dtype != 'category' and cleaned_data_dict[col][col].dtype != 'object' else None,
                    'showscale': False if cleaned_data_dict[col][col].dtype == 'category' or cleaned_data_dict[col][col].dtype == 'object' else True,
                    'line': {'width': 0.5, 'color': 'black'}
                },
            }, {
                'annotations[0].text': f"{col}: F={permanova_results[col]['test statistic']:.3f}, R<sup>2</sup>={permanova_results[col]['R²']:.3f}, p={permanova_results[col]['p-value']:.3f}" if col in permanova_results else ""
            }],
            label=col,
            method='update'
        ) for col in permanova_columns
    ]),
    direction='down',
    pad={'r': 10, 't': 10},
    showactive=True,
    x=1.1,
    xanchor='right',
    y=1.1,
    yanchor='top'
)

# Create base PCoA plot
initial_col = permanova_columns[0]
pcoa_fig = go.Figure(data=go.Scatter(
    x=cleaned_data_dict[initial_col]['PC1'],
    y=cleaned_data_dict[initial_col]['PC2'],
    mode='markers',
    marker=dict(
        color=cleaned_data_dict[initial_col][initial_col].astype('category').cat.codes if cleaned_data_dict[initial_col][initial_col].dtype == 'category' or cleaned_data_dict[initial_col][initial_col].dtype == 'object' else cleaned_data_dict[initial_col][initial_col],
        colorscale='Viridis' if cleaned_data_dict[initial_col][initial_col].dtype != 'category' and cleaned_data_dict[initial_col][initial_col].dtype != 'object' else None,
        showscale=False if cleaned_data_dict[initial_col][initial_col].dtype == 'category' or cleaned_data_dict[initial_col][initial_col].dtype == 'object' else True,
        line={'width': 0.5, 'color': 'black'}
    ),
    text=create_hover_text(cleaned_data_dict[initial_col]),
    hoverinfo='text'
))

pcoa_fig.update_layout(
    title=f'PCoA based on Bray-Curtis Distance',
    xaxis_title='PC1', yaxis_title='PC2',
    title_font_size=16, font_size=14, plot_bgcolor='white',
    updatemenus=[permanova_dropdown]
)

# Add initial PERMANOVA result as annotation if available
if initial_col in permanova_results:
    annotation = go.layout.Annotation(
        x=1, y=1, xref="paper", yref="paper", xanchor="right", yanchor="top",
        text=f"{initial_col} - PERMANOVA stats: F={permanova_results[initial_col]['test statistic']:.3f}, R²={permanova_results[initial_col]['R²']:.3f}, p={permanova_results[initial_col]['p-value']:.3f}",
        showarrow=False, font=dict(size=12), bgcolor="white", opacity=0.8
    )
    pcoa_fig.update_layout(annotations=[annotation])

# Save the plot as an HTML file
plot(pcoa_fig, filename='pcoa_plot_interactive.html', auto_open=False)
