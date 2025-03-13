import os 
import argparse
# Load gene count data and metadata
from plotly.subplots import make_subplots

parser = argparse.ArgumentParser(description='PCoA plot with PERMANOVA - metaFun')
parser.add_argument('--cpus', type=int, default=8, help='Number of CPU threads to use (default: 8)')
parser.add_argument('--metadata', type=str, required=True, help='Path to the metadata file')
from plotly.io import to_html

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
#metadata = pd.read_csv('selected_metadata.csv', index_col=0)
metadata = pd.read_csv(args.metadata, index_col=0)
custom_html = """
<style>
    .js-plotly-plot, .plot-container {
        width: 100%;
        height: 100%;
    }
    .scrollable-menu .Select-menu-outer {
        max-height: 300px !important;
        overflow-y: auto !important;
    }
</style>
<script>
    document.addEventListener('DOMContentLoaded', (event) => {
        const updateMenu = document.querySelector('.updatemenu-container');
        if (updateMenu) {
            updateMenu.classList.add('scrollable-menu');
            const dropdown = updateMenu.querySelector('.Select-menu-outer');
            if (dropdown) {
                dropdown.addEventListener('wheel', (e) => {
                    e.preventDefault();
                    dropdown.scrollTop += e.deltaY;
                }, { passive: false });
            }
        }
        
        // 그래프 크기 조정 함수
        function resizePlot() {
            const plot = document.querySelector('.js-plotly-plot');
            if (plot) {
                plot.style.height = window.innerHeight + 'px';
                Plotly.Plots.resize(plot);
            }
        }

        // 초기 크기 설정 및 창 크기 변경 시 대응
        resizePlot();
        window.addEventListener('resize', resizePlot);
    });
</script>
"""
config = {
    'responsive': True
}




# Calculate Bray-Curtis distance matrix and PCoA
bray_dm = beta_diversity("braycurtis", gene_counts.T)
bray_pcoa = pcoa(bray_dm)

# Prepare PCoA results for plotting
pcoa_results = pd.DataFrame(bray_pcoa.samples.iloc[:, :2], columns=['PC1', 'PC2'])
pcoa_results['sample_id'] = pcoa_results.index
merged_data = pcoa_results.merge(metadata, left_on='sample_id', right_index=True)

# Select columns for PERMANOVA
#permanova_columns = metadata.columns[np.r_[0, 1, 14:len(metadata.columns)]]
permanova_columns = metadata.columns[np.r_[0, :len(metadata.columns)]]

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
    if cleaned_data[col].nunique() > 1 and len(cleaned_data) >= 5:  # 최소 5개 샘플 필요

#    if cleaned_data[col].nunique() > 1:
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
            print(f"Error in PERMANOVA for {col}: {e}")
            print(f"Skipping {col} due to error: {e}")
    else:
        if cleaned_data[col].nunique() <= 1:
            print(f"Skipping {col} due to insufficient unique values after NA removal.")
        else:
            print(f"Skipping {col} due to insufficient samples (only {len(cleaned_data)} samples available, minimum 5 required).")        

# Create dropdown menu for PERMANOVA results
# def create_hover_text(df):
#     hover_text = df.apply(lambda row: '<br>'.join([f"{col}: {row[col]}" for col in df.columns]), axis=1)
#     return hover_text
def create_hover_text(df):
    def format_text(row):
        formatted_text = []
        for col in df.columns:
            text = str(row[col]).replace(';', '<br>')
            formatted_text.append(f"{col}: {text}")
        return '<br>'.join(formatted_text)

    hover_text = df.apply(format_text, axis=1)
    return hover_text

# permanova_dropdown = go.layout.Updatemenu(
#     buttons=list([
#         dict(
#             args=[{
#                 'marker': {
#                     'color': cleaned_data_dict[col][col].astype('category').cat.codes if cleaned_data_dict[col][col].dtype == 'category' or cleaned_data_dict[col][col].dtype == 'object' else cleaned_data_dict[col][col],
#                     'colorscale': 'Viridis' if cleaned_data_dict[col][col].dtype != 'category' and cleaned_data_dict[col][col].dtype != 'object' else None,
#                     'showscale': False if cleaned_data_dict[col][col].dtype == 'category' or cleaned_data_dict[col][col].dtype == 'object' else True,
#                     'name': col,
#                     'line': {'width': 0.5, 'color': 'black'}
#                 },
#             }, {
#                 'annotations[0].text': f"{col}: F={permanova_results[col]['test statistic']:.3f}, R<sup>2</sup>={permanova_results[col]['R²']:.3f}, p={permanova_results[col]['p-value']:.3f}" if col in permanova_results else ""
#             }],
#             label=col,
#             method='update'
#         ) for col in permanova_columns
#     ]),
#     direction='down',
#     pad={'r': 10, 't': 10},
#     showactive=True,
#     x=1.1,
#     xanchor='right',
#     y=1.1,
#     yanchor='top'
# )

#Create base PCoA plot
#initial_col = permanova_columns[0]
# pcoa_fig = go.Figure(data=go.Scatter(
#     x=cleaned_data_dict[initial_col]['PC1'],
#     y=cleaned_data_dict[initial_col]['PC2'],
#     mode='markers',
#     marker=dict(
#         color=cleaned_data_dict[initial_col][initial_col].astype('category').cat.codes if cleaned_data_dict[initial_col][initial_col].dtype == 'category' or cleaned_data_dict[initial_col][initial_col].dtype == 'object' else cleaned_data_dict[initial_col][initial_col],
#         colorscale='Viridis' if cleaned_data_dict[initial_col][initial_col].dtype != 'category' and cleaned_data_dict[initial_col][initial_col].dtype != 'object' else None,
#         showscale=False if cleaned_data_dict[initial_col][initial_col].dtype == 'category' or cleaned_data_dict[initial_col][initial_col].dtype == 'object' else True,
#         line={'width': 0.5, 'color': 'black'}
#     ),
#     text=create_hover_text(cleaned_data_dict[initial_col]),
#     hoverinfo='text'
# ))



# Create PCoA plot newly for metadata integration 
def create_pcoa_plot_for_column(fig, col):
    fig.data = [trace for trace in fig.data if not trace.name.startswith(f"{col}:") and trace.name != col]



    if cleaned_data_dict[col][col].dtype == 'bool' or cleaned_data_dict[col][col].dtype == 'category' or cleaned_data_dict[col][col].dtype == 'object':
        unique_categories = cleaned_data_dict[col][col].unique()
        for category in unique_categories:
            category_text = category if isinstance(category, str) else str(category)
            category_text = category_text.replace(';', '<br>')
            category_data = cleaned_data_dict[col][cleaned_data_dict[col][col] == category]
            fig.add_trace(go.Scatter(
                x=category_data['PC1'],
                y=category_data['PC2'],
                mode='markers',
                name=f"{col}:<br>{category_text}",  # Updated to add line breaks in the legend
                marker=dict(
                    size=7, opacity=0.7,
                    line=dict(width=0.5, color='black')
                ),
                text=create_hover_text(category_data),
                hoverinfo='text',
                showlegend=True  # Ensure legend is shown
            ))
    else:
        fig.add_trace(go.Scatter(
            x=cleaned_data_dict[col]['PC1'],
            y=cleaned_data_dict[col]['PC2'],
            mode='markers',
            name=col,
            marker=dict(
                color=cleaned_data_dict[col][col],
                colorscale='Viridis',
                showscale=True,
                size=7, opacity=0.7,
                colorbar=dict(y=0.45,
                 #   title=col,
                    titleside='top'  # Position the title above the color bar
                ),                
                line={'width': 0.5, 'color': 'black'}
            ),
            text=create_hover_text(cleaned_data_dict[col]),
            hoverinfo='text',
            showlegend=True  # Ensure legend is shown
        ))

pcoa_fig = go.Figure()




pcoa_fig.add_layout_image(
    dict(
        xref="paper", yref="paper",
        x=0.5, y=0.5,
        sizex=1, sizey=1,
        xanchor="center", yanchor="middle"
    )
)



# Ensure initial column is the first column
last_col = permanova_columns[-1]
initial_col =  permanova_columns[0]

#create_pcoa_plot_for_column(pcoa_fig, initial_col)

# for col in permanova_columns:
#     create_pcoa_plot_for_column(pcoa_fig, col)
#     pcoa_fig.data[-1].visible = False  # Only initial column traces visible

#pcoa_fig.data[0].visible = False 


# for col in permanova_columns:
#     create_pcoa_plot_for_column(pcoa_fig, col)
#     if col == last_col:
#         for trace in pcoa_fig.data[-len(cleaned_data_dict[col][col].unique()):]:
#             trace.visible = True
#     else:
#         for trace in pcoa_fig.data[-len(cleaned_data_dict[col][col].unique()):]:
#             trace.visible = False

for col in permanova_columns:
    create_pcoa_plot_for_column(pcoa_fig, col)
    # 초기 컬럼(Completeness)을 제외한 모든 트레이스를 숨김
    if col != initial_col:
        for trace in pcoa_fig.data[-len(cleaned_data_dict[col][col].unique()):]:
            trace.visible = False

# 초기 컬럼(Completeness)의 트레이스만 표시
for trace in pcoa_fig.data:
    if trace.name.startswith(f"{initial_col}:") or trace.name == initial_col:
        trace.visible = True



genome_data_columns = permanova_columns[:21]
metadata_columns = permanova_columns[21:]

# def create_dropdown_menu(columns, name, x_position, is_genome_data):
#     buttons = []
#     for col in columns:
#         buttons.append(
#             dict(
#                 args=[{
#                     "visible": [
#                         (trace.name.startswith(f"{col}:") or trace.name == col)
#                         for trace in pcoa_fig.data
#                     ] + [False] * len(pcoa_fig.data),
#                     "showlegend": [
#                         (trace.name.startswith(f"{col}:") or trace.name == col)
#                         for trace in pcoa_fig.data
#                     ] + [False] * len(pcoa_fig.data),
#                     "updatemenus[0].active": 0 if not is_genome_data else None,
#                     "updatemenus[1].active": 0 if is_genome_data else None
#                 }, {
#                     "annotations[1].text": f"<span style='color:blue;'>{col} : PERMANOVA stats: F={permanova_results[col]['test statistic']:.3f}, R²={permanova_results[col]['R²']:.3f}, p={permanova_results[col]['p-value']:.3f}</span>" if col in permanova_results else ""
#                 }],
#                 label=col,
#                 method="update"
#             )
#         )
    
#     # Add a "None" option at the beginning of the dropdown
#     buttons.insert(0, dict(
#         args=[{
#             "visible": [False] * len(pcoa_fig.data),
#             "showlegend": [False] * len(pcoa_fig.data),
#             "updatemenus[0].active": 0 if is_genome_data else None,
#             "updatemenus[1].active": 0 if not is_genome_data else None
#         }, {
#             "annotations[1].text": ""
#         }],
#         label="None",
#         method="update"
#     ))
    
#     return go.layout.Updatemenu(
#         buttons=buttons,
#         direction='down',
#         pad={'r': 10, 't': 10},
#         showactive=True,
#         x=x_position,
#         xanchor='left',
#         y=1.07,
#         yanchor='top',
#         name=name,
#         active=0  # Set the default selection to "None"
#     )
buttons = []
for col in genome_data_columns:
    buttons.append(dict(
        label=f"{col}",
        method="update",
        args=[{
            "visible": [trace.name.startswith(f"{col}:") or trace.name == col for trace in pcoa_fig.data],
            "showlegend": [trace.name.startswith(f"{col}:") or trace.name == col for trace in pcoa_fig.data]
        }, {
            "annotations[1].text": f"<span style='color:blue;'>{col} : PERMANOVA stats: F={permanova_results[col]['test statistic']:.3f}, R²={permanova_results[col]['R²']:.3f}, p={permanova_results[col]['p-value']:.3f}</span>" if col in permanova_results else ""
        }]
    ))

# Metadata 버튼
for col in metadata_columns:
    buttons.append(dict(
        label=f"{col}",
        method="update",
        args=[{
            "visible": [trace.name.startswith(f"{col}:") or trace.name == col for trace in pcoa_fig.data],
            "showlegend": [trace.name.startswith(f"{col}:") or trace.name == col for trace in pcoa_fig.data]
        }, {
            "annotations[1].text": f"<span style='color:blue;'>{col} : PERMANOVA stats: F={permanova_results[col]['test statistic']:.3f}, R²={permanova_results[col]['R²']:.3f}, p={permanova_results[col]['p-value']:.3f}</span>" if col in permanova_results else ""
        }]
    ))


# genome_data_dropdown = create_dropdown_menu(genome_data_columns, "Genome Data", 0.5, True)
# metadata_dropdown = create_dropdown_menu(metadata_columns, "Metadata", 0.92, False)

pcoa_fig.data[0].visible = True
for trace in pcoa_fig.data[1:]:
    trace.visible = False
# newly added for button 
########################


link_text= 'Click here for more information about metaFun'
link_url = 'https://metafun-doc.readthedocs.io/en/latest/workflows/COMPARATIVE_ANNOTATION.html'

link_annotation = go.layout.Annotation(
    xref="paper", yref="paper", xanchor="left", yanchor="top", x=0, y=1.05,
    text=f'<a href="{link_url}" style="color:blue;">{link_text}</a>',
    showarrow=False, font=dict(size=14), bgcolor="white", opacity=0.8
)

annotation = None  

# Add initial PERMANOVA result as annotation if available
if initial_col in permanova_results:
    annotation = go.layout.Annotation(
        x=0.1, y=0.1, xref="paper", yref="paper", xanchor="left", yanchor="top",
        text=f"<span style='color:blue;'>{initial_col} : PERMANOVA stats: F={permanova_results[initial_col]['test statistic']:.3f}, R²={permanova_results[initial_col]['R²']:.3f}, p={permanova_results[initial_col]['p-value']:.3f}</span>",
        showarrow=False, font=dict(size=12), bgcolor="white", opacity=0.8
    )

else:
    annotation = go.layout.Annotation(
        x=0.1, y=0.1, xref="paper", yref="paper", xanchor="left", yanchor="top",
        text=f"<span style='color:red;'>{initial_col} : PERMANOVA skipped - insufficient data</span>",
        showarrow=False, font=dict(size=12), bgcolor="white", opacity=0.8
    )

dropdown_annotations = [
    dict(
        text="Genome Data",
        x=0.75,
        xref="paper",
        y=1.12,
        yref="paper",
        showarrow=False,
        align="center"
    ),
    dict(
        text="Metadata",
        x=1.05,
        xref="paper",
        y=1.12,
        yref="paper",
        showarrow=False,
        align="center"
    )
]

# annotations_list = [link_annotation]
# if permanova_annotation:
#     annotations_list.append(permanova_annotation)

pcoa_fig.update_layout(
    title=f'PCoA based on Bray-Curtis Distance',
    xaxis_title='PC1', yaxis_title='PC2',
    title_font_size=17, font_size=14, plot_bgcolor='white',
    #updatemenus=[permanova_dropdown],
    updatemenus=[dict(
        buttons=buttons,
        direction="down",
        pad={"r": 10, "t": 10},
        showactive=True,
        x=0.9,
        xanchor="left",
        y=1.15,
        yanchor="top",
    )],
    annotations=[link_annotation, annotation],# +  dropdown_annotations,
    #width=900,
    #height=600,
    #
    margin=dict(t=100),
    xaxis=dict(
        title='PC1',
        showgrid=True, 
        gridcolor='lightgray',
        gridwidth=0.5,  
        griddash='dash' ,
        zeroline=True,  
        zerolinecolor='lightgray',
        zerolinewidth=0.5 
    ),
    yaxis=dict(
        title='PC2',
        showgrid=True,  
        gridcolor='lightgray',
        gridwidth=0.5,  
        griddash='dash' ,
        zeroline=True,  
        zerolinecolor='lightgray',
        zerolinewidth=0.5 
    ),    
)

initial_col = genome_data_columns[0]
for trace in pcoa_fig.data:
    trace.visible = trace.name.startswith(f"{initial_col}:") or trace.name == initial_col



# pcoa_fig.update_layout(
#     margin=dict(t=100),
#     updatemenus=[{
#         "buttons": updatemenu_buttons,
#         "direction": 'down',
#         "pad": {'r': 10, 't': 10},
#         "showactive": True,
#         "x": 1.0,
#         "xanchor": 'left',
#         "y": 1.1,
#         "yanchor": 'top',
#         "active": len(permanova_columns) - 1 # last column 
#     }]
# )


#    pcoa_fig.update_layout(annotations=[annotation])

# Save the plot as an HTML file
#plot(pcoa_fig, filename='pcoa_plot_interactive.html', auto_open=False)

plot_html = to_html(pcoa_fig, include_plotlyjs=True, full_html=False)
full_html = f"""
<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>PCoA Plot</title>
    {custom_html}
</head>
<body>
    {plot_html}
</body>
</html>
"""


with open('pcoa_plot_interactive.html', 'w', encoding='utf-8') as f:
    f.write(full_html)