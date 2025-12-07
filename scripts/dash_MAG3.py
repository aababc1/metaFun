# good sunburst 
import dash
from dash import dcc, html, dash_table
from dash.dependencies import Input, Output, State
from dash import callback_context
import plotly.express as px
import pandas as pd
import plotly.graph_objects as go
from dash.exceptions import PreventUpdate
import os 
import numpy as np
import argparse 
import dash_daq as daq 
import dash_bootstrap_components as dbc

parser = argparse.ArgumentParser(description='metaFun: genome selector for COMPARATIVE_ANNOTATION')
parser = argparse.ArgumentParser(description='refer to the documentation at https://metafun-doc-v01.readthedocs.io/en/latest/')
parser.add_argument('-i', '--input', help='Input CSV file', required=True)
parser.add_argument('-p', '--port', help='Port to run the server on', type=int, default=8050)
args = parser.parse_args()


#app = dash.Dash(__name__, suppress_callback_exceptions=True)
app = dash.Dash(__name__, external_stylesheets=[dbc.themes.BOOTSTRAP], suppress_callback_exceptions=True)

# Read the data from the CSV file
df = pd.read_csv(args.input)
df_raw = pd.read_csv(args.input)


#df = df.fillna('NA')
#df = df.apply(lambda col: col.fillna(-999) if col.dtype.kind in 'biufc' else col.fillna('Missing'))
df = df.fillna('None')
"""
def prepare_data(df):
    for col in df.columns:
        if pd.api.types.is_numeric_dtype(df[col]):
            df[col] = df[col].fillna(-1)  # Or another appropriate marker
        else:
            df[col] = df[col].fillna('None')
    return df

df = prepare_data(df)
"""
def create_scatter_plot(df, x_col, y_col, color_col):
    # Mask to identify rows with placeholder values
    mask = df[color_col] == 'None'

    # Scatter plot for regular data
    scatter_fig = px.scatter(
        df[~mask],
        x=x_col,
        y=y_col,
        color=color_col,
        labels={"color": color_col},
        title='Completeness vs Contamination',
        color_continuous_scale=px.colors.sequential.Viridis
    )

    # Add traces for placeholder values, if any exist
    if mask.any():
        scatter_fig.add_trace(go.Scatter(
            x=df[mask][x_col],
            y=df[mask][y_col],
            mode='markers',
            marker=dict(color='red', size=10, symbol='x'),
            name='Missing'
        ))

    scatter_fig.update_layout(plot_bgcolor='white')
    return scatter_fig


# reformat classificaiton  for readibility 
def format_classification(text, break_after=3):
    # Detect the delimiter by checking the occurrence
    delimiter = ';' if ';' in text else '/'
    
    parts = text.split(delimiter)
    if len(parts) > break_after:
        # Insert a line break after the third part
        parts[break_after] = "<br>" + parts[break_after]
    return delimiter.join(parts)
    
df['formatted_classification'] = df['classification'].apply(format_classification)

# check whether the column is numeric
def is_numeric(series):
    """ Check if a pandas series is numeric. """
    return pd.api.types.is_numeric_dtype(series)

# Define taxonomic ranks
tax_ranks = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species']

# Split the classification column into separate columns for each level
df[tax_ranks] = df['classification'].str.split(';', expand=True)

metadata_columns = [col for col in df.columns if col not in tax_ranks + ['classification', 'Genome']]
dropdown_options = [{'label': col, 'value': col} for col in metadata_columns if col != 'formatted_classification']





def create_sunburst(df, selected_metadata):
    if df is None or df.empty or selected_metadata not in df.columns:
        fig = go.Figure()
        fig.add_annotation(text="data is not valid", 
                          xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        fig.update_layout(height=700, width=700)
        return fig
    
    try:
        #  work_df is a copy of df
        work_df = df.copy()
        
        # taxonomy data preprocessing - None values are converted to empty strings
        for rank in tax_ranks:
            work_df[rank] = work_df[rank].fillna('').replace('None', '')
        
        # Processing of selected metadata
        if pd.api.types.is_numeric_dtype(df[selected_metadata].dropna()):
            # Numeric type 
            # 1. None values are converted to NaN
            work_df[selected_metadata] = pd.to_numeric(work_df[selected_metadata], errors='coerce')
            
            # 2. Remove rows with NaN values
            valid_data = work_df.dropna(subset=[selected_metadata])
            
            # 3. Remove rows with empty taxonomy paths
            valid_data = valid_data[valid_data[tax_ranks].ne('').any(axis=1)]
            
            if valid_data.empty:
                fig = go.Figure()
                fig.add_annotation(text="No valid data after filtering NA values", 
                                  xref="paper", yref="paper",
                                  x=0.5, y=0.5, showarrow=False)
                fig.update_layout(height=700, width=700)
                return fig
            
            # 4. Binning for numeric data - This is the key point
            # Binning is used to categorize numeric data
            n_bins = min(10, len(valid_data[selected_metadata].unique()))
            if n_bins <= 1:
                n_bins = 2  # Use at least 2 bins
            
            # Create binned categories
            valid_data['binned'] = pd.cut(
                valid_data[selected_metadata], 
                bins=n_bins,
                labels=[f"{selected_metadata}: {i}" for i in range(n_bins)]
            )
            
            # Set path
            path = [col for col in tax_ranks if valid_data[col].nunique() > 1] + ['binned']
            
            # Create sunburst plot when color is numeric data
            fig = px.sunburst(
                valid_data,
                path=path,
                color=selected_metadata,
                color_continuous_scale=px.colors.sequential.RdBu
            )
        else:
            # String/categorical data processing
            # Remove rows with None values
            valid_data = work_df[work_df[selected_metadata].notna() & 
                               (work_df[selected_metadata] != 'None') & 
                               (work_df[selected_metadata] != '')]
            
            # Remove rows with empty taxonomy paths
            valid_data = valid_data[valid_data[tax_ranks].ne('').any(axis=1)]
            
            if valid_data.empty:
                fig = go.Figure()
                fig.add_annotation(text="No valid data after filtering NA values", 
                                  xref="paper", yref="paper",
                                  x=0.5, y=0.5, showarrow=False)
                fig.update_layout(height=700, width=700)
                return fig
            
            # Only include paths that are actually used (columns with only one value are not meaningful)
            path = [col for col in tax_ranks if valid_data[col].nunique() > 1] + [selected_metadata]
            
            fig = px.sunburst(
                valid_data,
                path=path,
                color=selected_metadata
            )
        
        # Layout setting
        fig.update_layout(
            height=800, 
            width=800,
            margin=dict(t=30, l=0, r=0, b=0)
        )
        
        return fig
        
    except Exception as e:
        # When an error occurs, print the error and provide feedback
        print(f"Sunburst plot creation error: {str(e)}")
        fig = go.Figure()
        fig.add_annotation(text=f"Sunburst plot creation error: {str(e)}", 
                          xref="paper", yref="paper",
                          x=0.5, y=0.5, showarrow=False)
        fig.update_layout(height=700, width=700)
        return fig


section_style = {
    'backgroundColor': '#f0f0f0',
    'padding': '20px',
    'margin': '10px 0',
    'borderRadius': '5px',
    'boxShadow': '0 2px 4px rgba(0,0,0,0.1)'
}
common_style = {
    'border': '2px solid #ddd',
    'borderRadius': '5px',
    'padding': '20px',
    'marginBottom': '20px',
    'backgroundColor': 'white',
}
### app.layout region 
### app.layout region 
app.layout = html.Div([
    html.H1('metaFun : genome selector for COMPARATIVE_ANNOTATION', style={'textAlign': 'center'}),

    html.Div([
        html.Div([
            html.H2('Quality and number of all genomes', style={'textAlign': 'left', 'marginBottom': '20px'}),
            html.P([
                "This section shows the quality metrics of all genomes and their distribution. ",
                "The left plot displays genome quality as Completeness vs Contamination. ",
                "The right plot shows the distribution of selected metadata values.",
                "You can select your metadata or taxonomic rank to color and identify distribution of your genomes."
            ], style={
                'color': '#0066cc', 
                'fontSize': '14px', 
                'marginBottom': '15px',
                'fontStyle': 'italic'
            }),
            html.Div([
                html.Div([
                    html.Label('Select your metadata', style={'fontSize': 15, 'fontWeight': 'bold'}),
                    dcc.Dropdown(
                        id='color-dropdown',
                        options=dropdown_options,
                        value=metadata_columns[0] if metadata_columns else None,
                        placeholder="Select metadata",
                        style={'width': '200px'}
                    ),
                    dbc.Tooltip(
                        "Select metadata to color and analyze your genomes",
                        target="color-dropdown",
                        style={"background-color": "#17a2b8", "color": "white"}
                    ),
                ], style={'display': 'inline-block', 'marginRight': '20px', 'verticalAlign': 'top'}),
                html.Div([
                    html.Label('Select taxonomy rank', style={'fontSize': 15, 'fontWeight': 'bold'}),
                    dcc.Dropdown(
                        id='tax-rank-dropdown',
                        options=[{'label': rank, 'value': rank} for rank in tax_ranks],
                        value=tax_ranks[0],
                        style={'width': '200px'}
                    ),
                    dbc.Tooltip(
                        "Select taxonomic rank to filter genomes (domain to species)",
                        target="tax-rank-dropdown",
                        style={"background-color": "#17a2b8", "color": "white"}
                    ),
                ], style={'display': 'inline-block', 'marginRight': '20px', 'verticalAlign': 'top'}),
                html.Div([
                    html.Label('Select taxonomy values', style={'fontSize': 15, 'fontWeight': 'bold'}),
                    dcc.Dropdown(
                        id='tax-value-dropdown',
                        multi=True,
                        style={'width': '200px'}
                    ),
                    dbc.Tooltip(
                        "Select specific taxonomic groups to analyze",
                        target="tax-value-dropdown",
                        style={"background-color": "#17a2b8", "color": "white"}
                    ),
                ], style={'display': 'inline-block', 'verticalAlign': 'top'}),
            ], style={'marginBottom': '20px', 'textAlign': 'right'}),
        ], style={'backgroundColor': '#f0f0f0', 'padding': '20px', 'borderRadius': '5px', 'marginBottom': '20px', 'width': '100%', 'display': 'inline-block'}),

        html.Div([
            html.H3('Scatter Plot of Genome Quality', style={'textAlign': 'center'}), 
            dcc.Graph(
                id='overall-scatter-plot', 
                selectedData={'points': []},    
                config={
                    'displayModeBar': True,
                    'modeBarButtonsToAdd': ['select2d', 'lasso2d'],
                    'scrollZoom': True
                },
                style={'height': '600px'}
            ),
            dbc.Tooltip(
                "Scatter plot showing genome quality metrics. Use lasso or box select to analyze specific genomes",
                target="overall-scatter-plot",
                placement="top",
                style={"background-color": "#17a2b8", "color": "white"}
            )
        ], style={'width': '65%', 'display': 'inline-block', 'overflow': 'visible'}),
        
        html.Div([
            html.H3('Plot of Selected Metadata', style={'textAlign': 'center'}), 
            dcc.Graph(
                id='overall-distribution-plot', 
                style={'height':'600px'}
            ),
            dbc.Tooltip(
                "Distribution of selected metadata across genomes",
                target="overall-distribution-plot",
                placement="top",
                style={"background-color": "#17a2b8", "color": "white"}
            )
        ], style={'width': '35%', 'display': 'inline-block', 'vertical-align': 'top', 'overflow': 'visible'}),
    ], style={**common_style, 'overflow': 'visible'}),    

    html.Div([
        html.Div([
            html.H2('Select genomes that you are interested in', style={'textAlign': 'left'}),


        html.P([
            "Clicking on a section of the Sunburst Plot will filter genomes by taxonomy. "
            "The quality information of selected genomes will be visualized below. "
            "You can download only these selected genomes using the download buttons.",

        ], style={
                        'color': '#0066cc', 
                        'fontSize': '14px', 
                        'marginBottom': '15px',
                        'fontStyle': 'italic'
                    }),
        html.P([
            "Result csv file is utilized in COMPARATIVE_ANNOTATION module."
        ], style={
                        'color': '#0066cc', 
                        'fontSize': '14px', 
                        'marginBottom': '15px',
                        'fontStyle': 'italic'
                    }),
            html.Label('Select your metadata for subsetting genomes', 
                      style={'fontSize': 15, 'fontWeight': 'bold', 'marginBottom': '10px'}),
            html.Div([
                dcc.Dropdown(
                    id='metadata-dropdown',
                    options=dropdown_options,
                    value=metadata_columns[0] if metadata_columns else None,
                    style={'width': '300px'}
                ),
                dbc.Tooltip(
                    "Select metadata to analyze in the sunburst plot",
                    target="metadata-dropdown",
                    placement='bottom',
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
                html.Div(id='dynamic-filtering')
            ],style={ 'width':'100%','margin': '10px auto','textAlign': 'left'}),
        ], style=section_style),

        html.Div([
            html.H3('Sunburst Plot Displaying Taxonomic Ranks and Selected Metadata', style={'textAlign': 'center', 'marginBottom': '10px'}),            
            html.Div([
                html.Button('Reset to Top Level', 
                          id='reset-sunburst-button', 
                          n_clicks=0),
                dbc.Tooltip(
                    "Reset sunburst plot to show all data",
                    target="reset-sunburst-button",
                    placement='bottom',
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
                html.Div(id='sunburst-level-display', 
                        children='You are seeing now: Entire Dataset')
            ], style={'margin-bottom': '10px'}),
            html.Div([
                dcc.Graph(id='sunburst-plot'),
                dbc.Tooltip(
                    "Interactive visualization of genome distribution. Click segments to zoom in",
                    target="sunburst-plot",
                    placement="top",
                    style={"background-color": "#17a2b8", "color": "white"}
                )
            ], style={'display': 'flex', 'justify-content': 'center', 'align-items': 'center', 'height': '100%'}),
        ], style={'overflow': 'visible'}),

        html.Div([
            # left scatter plot
            html.Div([
                html.H3('Genome Quality Scatter Plot of Selected Taxon', 
                       style={'textAlign': 'center', 'marginBottom': '10px'}),
                dcc.Graph(
                    id='filtered-scatter-plot', 
                    style={'height': '600px'}
                ),
                dbc.Tooltip(
                    "Filtered scatter plot based on selected criteria",
                    target="filtered-scatter-plot",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
            ], style={'width': '65%', 'display': 'inline-block'}),
            
            # right distribution plot
            html.Div([
                html.H3('Metadata Distribution Plot of Selected Taxon', 
                       style={'textAlign': 'center', 'marginBottom': '10px'}),
                dcc.Graph(
                    id='filtered-distribution-plot', 
                    style={'height': '600px'}
                ),
                dbc.Tooltip(
                    "Distribution of filtered data",
                    target="filtered-distribution-plot",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
            ], style={'width': '35%', 'display': 'inline-block'}),
        ], style={'width': '100%', 'display': 'flex', 'overflow': 'visible'}),

        # added start 
        dbc.Modal(
            [
                dbc.ModalHeader("File Already Exists"),
                dbc.ModalBody("A file with this name already exists. Do you want to overwrite it?"),
                dbc.ModalFooter([
                    dbc.Button("Cancel", id="cancel-overwrite", className="ml-auto"),
                    dbc.Button("Overwrite", id="confirm-overwrite", color="danger"),
                ]),
            ],
            id="overwrite-modal",
            is_open=False,
            style={
                'zIndex': 9999,  # higher z-index value to display above other elements
                'backgroundColor': 'rgba(0, 0, 0, 0.5)'  # darken the background
            },
            centered=True  # display in the center

        ),
        html.Div(id='temp-storage', style={'display': 'none'}),
        #added end

        html.Div([
            html.H3('Download Metadata Table for COMPARATIVE_ANNOTATION', 
                   style={'textAlign': 'left', 'marginBottom': '10px'}),
                   
            html.Div([
                dcc.Input(
                    id='search-input',
                    type='text',
                    placeholder='Search table...',
                    style={'marginRight': '10px', 'width': '200px'}
                ),
                dbc.Tooltip(
                    "Search for specific terms in the table",
                    target="search-input",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
                dcc.Input(
                    id='filename-input',
                    type='text',
                    placeholder='Enter filename for CSV',
                    value='genome_selector_result',  
                    
                    style={'marginRight': '10px', 'width': '200px'}
                ),
                dbc.Tooltip(
                    "Enter filename to save your selected data",
                    target="filename-input",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
                html.Button('Save to Local', id='save-local-button', n_clicks=0),
                dbc.Tooltip(
                    "Download selected data as CSV file",
                    target="save-local-button",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
                html.Button('Save to Server', id='save-server-button', n_clicks=0),
                dbc.Tooltip(
                    "Save selected data to server",
                    target="save-server-button",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
            ], style={'marginBottom': '10px', 'marginTop': '20px'}),
        



            html.H4('Selected Part of Metadata Table', 
                style={'textAlign': 'left', 'marginTop': '20px', 'marginBottom': '10px'}),        
            html.Div([
                dcc.Dropdown(
                    id='column-selector',
                    options=[{'label': col, 'value': col} for col in df.columns if col != 'formatted_classification'],
                    value=[col for col in df.columns[:10]],
                    multi=True,
                    placeholder="Select columns to display"
                ),
                dbc.Tooltip(
                    "Choose which columns to display in the table",
                    target="column-selector",
                    style={"background-color": "#17a2b8", "color": "white"}
                ),
            ], style={'width': '50%', 'margin': '10px auto'}),

            dash_table.DataTable(
                id='data-table',
                columns=[{"name": i, "id": i} for i in df.columns[:10] if i != 'formatted_classification'],  # Adjust number of columns displayed        
                data=df.head(20).to_dict('records'),
                page_size=20,
                style_table={'overflowX': 'scroll'},
                css=[{ 
                    'selector': '.dash-spreadsheet-page-selector',
                    'rule': 'display: flex; justify-content: flex-start; margin-left: 10px;'

                    }],  # Enable horizontal scroll
                filter_action='custom',
                filter_query='',
                sort_action='custom',
                sort_mode='multi',
                sort_by=[]
            ),
            dcc.Download(id="download-dataframe-csv"),
            html.Div(id='save-status') 
        ], style=common_style)
    ])
])
### app.layout region 
### app.layout region 
###################

@app.callback(
    Output('dynamic-filtering', 'children'),
    [Input('metadata-dropdown', 'value')]
)
def display_dynamic_filter(selected_metadata):
    # Do not display the filtering component and return an empty div
    return html.Div()


@app.callback(
    Output('tax-value-dropdown', 'options'),
    [Input('tax-rank-dropdown', 'value')]
)
def update_tax_value_options(selected_tax_rank):
    if not selected_tax_rank:
        return []
    unique_values = df[selected_tax_rank].unique()
    return [{'label': val, 'value': val} for val in unique_values if val != 'None']

@app.callback(
    [Output('overall-scatter-plot', 'figure'),
     Output('overall-distribution-plot', 'figure')],
    [Input('color-dropdown', 'value'),
     Input('tax-rank-dropdown', 'value'),
     Input('tax-value-dropdown', 'value'),
     Input('overall-scatter-plot', 'selectedData')]
)
def update_overall_plots(color_var, selected_tax_rank, selected_tax_values, selected_data):
    # Create a copy of the dataframe for plot manipulation
    plot_df = df.copy()
    
    # Convert 'None' placeholders to -999 for numeric data
    plot_df[color_var] = plot_df[color_var].replace('None', -999)

    if selected_tax_rank and selected_tax_values:
        plot_df = plot_df[plot_df[selected_tax_rank].isin(selected_tax_values)]
        
    # Mask for identifying placeholders in the dataset
    placeholder_mask = plot_df[color_var] == -999

    if selected_data and selected_data['points']:
        selected_points = selected_data['points']
        selected_indices = [point['pointIndex'] for point in selected_points]
        filtered_data = plot_df.iloc[selected_indices]
    else:
        filtered_data = plot_df

    # Create scatter plot with go.Figure 
    scatter_fig = go.Figure()

    # Check if the color variable contains numeric data
    is_numeric_data = pd.api.types.is_numeric_dtype(filtered_data[~placeholder_mask][color_var]) 

    # For non-numeric data, use a discrete colormap
    if not is_numeric_data:
        # Get unique categories and create a categorical color mapping
        if not filtered_data[~placeholder_mask].empty:
            unique_categories = filtered_data[~placeholder_mask][color_var].unique()
            # Use px.colors functions to generate color map
            category_colors = px.colors.qualitative.Plotly[:len(unique_categories)]
            color_map = dict(zip(unique_categories, category_colors))
            
            # Add a trace for each category
            for category in unique_categories:
                category_data = filtered_data[~placeholder_mask][filtered_data[~placeholder_mask][color_var] == category]
                scatter_fig.add_trace(go.Scatter(
                    x=category_data['Completeness'],
                    y=category_data['Contamination'],
                    mode='markers',
                    marker=dict(color=color_map[category]),
                    name=str(category),
                    text=category_data['formatted_classification'],
                    hovertemplate='Completeness: %{x}<br>Contamination: %{y}<br>Classification: %{text}<br>' + 
                                 f'{color_var}: {category}<extra></extra>'
                ))
    else:
        # For numeric data, use the original approach with a continuous colorscale
        if not filtered_data[~placeholder_mask].empty:
            scatter_fig.add_trace(go.Scatter(
                x=filtered_data[~placeholder_mask]['Completeness'],
                y=filtered_data[~placeholder_mask]['Contamination'],
                mode='markers',
                marker=dict(
                    color=filtered_data[~placeholder_mask][color_var].astype(float),
                    colorscale='Viridis',
                    showscale=True
                ),
                text=filtered_data[~placeholder_mask]['formatted_classification'],
                hovertemplate='Completeness: %{x}<br>Contamination: %{y}<br>Classification: %{text}<br>' + 
                             f'{color_var}: ' + '%{marker.color}<extra></extra>'
            ))

    # Add trace for placeholder values
    if placeholder_mask.any():
        scatter_fig.add_trace(go.Scatter(
            x=filtered_data[placeholder_mask]['Completeness'],
            y=filtered_data[placeholder_mask]['Contamination'],
            mode='markers',
            marker=dict(
                color='red',
                size=10,
                symbol='x',
                line=dict(color='Black', width=1)
            ),
            name='Not available',
            text=filtered_data[placeholder_mask]['formatted_classification'],
            hovertemplate='Completeness: %{x}<br>Contamination: %{y}<br>Classification: %{text}<extra></extra>'
        ))

    scatter_fig.update_layout(
        title='Overall: Completeness vs Contamination',
        plot_bgcolor='white',
        font={'family': 'Arial', 'size': 14},
        dragmode='select',
        clickmode='event+select',
        selectdirection='any',
        legend=dict(
            title='',
            orientation='h',
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        xaxis_title="Completeness",
        yaxis_title="Contamination"
    )

    # Distribution plot
    if selected_data and selected_data['points']:
        selected_indices = [point['pointIndex'] for point in selected_data['points']]
        df_for_hist = plot_df.iloc[selected_indices]
    else:
        df_for_hist = plot_df

    df_for_hist[color_var] = df_for_hist[color_var].replace(-999, 'Not available')
    
    if is_numeric_data:
        dist_fig = px.histogram(
            df_for_hist[df_for_hist[color_var] != 'Not available'],
            x=color_var,
            title=f'Distribution of {color_var} for Selected Points',
            nbins=30
        )
    else:
        value_counts = df_for_hist[color_var].value_counts()
        dist_fig = px.bar(
            x=value_counts.index,
            y=value_counts.values,
            title=f'Distribution of {color_var} for Selected Points'
        )
        dist_fig.update_xaxes(title_text=color_var)
        dist_fig.update_yaxes(title_text='Number of genomes')

    not_available_count = (df_for_hist[color_var] == 'Not available').sum()
    if not_available_count > 0:
        dist_fig.add_annotation(
            x=0.5, y=1.05,
            xref='paper', yref='paper',
            text=f'Number of genomes without metadata: {not_available_count}',
            showarrow=False,
            yshift=10,
            font=dict(color='red', size=12)
        )

    dist_fig.update_layout(plot_bgcolor='white', font={'family': 'Arial', 'size': 14})

    return scatter_fig, dist_fig
# @app.callback(
#     [Output('overall-scatter-plot', 'figure'),
#      Output('overall-distribution-plot', 'figure')],
#     [Input('color-dropdown', 'value'),
#      Input('tax-rank-dropdown', 'value'),
#      Input('tax-value-dropdown', 'value'),
#      Input('overall-scatter-plot', 'selectedData')]
# )

# def update_overall_plots(color_var, selected_tax_rank, selected_tax_values, selected_data):
#     # Create a copy of the dataframe for plot manipulation
#     plot_df = df.copy()
    
#     # Convert 'None' placeholders to -999 only for the plotting data
#     plot_df[color_var] = plot_df[color_var].replace('None', -999)

#     if selected_tax_rank and selected_tax_values:
#         plot_df = plot_df[plot_df[selected_tax_rank].isin(selected_tax_values)]
        
#     # Mask for identifying placeholders in the dataset
#     placeholder_mask = plot_df[color_var] == -999

#     if selected_data and selected_data['points']:
#         selected_points = selected_data['points']
#         selected_indices = [point['pointIndex'] for point in selected_points]
#         filtered_data = plot_df.iloc[selected_indices]
#     else:
#         filtered_data = plot_df

#     # Create scatter plot with go.Figure instead of px.scatter
#     scatter_fig = go.Figure()

#     # Add trace for non-placeholder data
#     if not filtered_data[~placeholder_mask].empty:
#         scatter_fig.add_trace(go.Scatter(
#             x=filtered_data[~placeholder_mask]['Completeness'],
#             y=filtered_data[~placeholder_mask]['Contamination'],
#             mode='markers',
#             marker=dict(
#                 color=filtered_data[~placeholder_mask][color_var],
#                 colorscale='Viridis',
#                 showscale=True
#             ),
#             text=filtered_data[~placeholder_mask]['formatted_classification'],
#             hovertemplate='Completeness: %{x}<br>Contamination: %{y}<br>Classification: %{text}<br>' + f'{color_var}: ' + '%{marker.color}<extra></extra>'
#         ))

#     # Add trace for placeholder values
#     if placeholder_mask.any():
#         scatter_fig.add_trace(go.Scatter(
#             x=filtered_data[placeholder_mask]['Completeness'],
#             y=filtered_data[placeholder_mask]['Contamination'],
#             mode='markers',
#             marker=dict(
#                 color='red',
#                 size=10,
#                 symbol='x',
#                 line=dict(color='Black', width=1)
#             ),
#             name='Not available',
#             text=filtered_data[placeholder_mask]['formatted_classification'],
#             hovertemplate='Completeness: %{x}<br>Contamination: %{y}<br>Classification: %{text}<extra></extra>'
#         ))

#     scatter_fig.update_layout(
#         title='Overall: Completeness vs Contamination',
#         plot_bgcolor='white',
#         font={'family': 'Arial', 'size': 14},
#         dragmode='select',
#         clickmode='event+select',
#         selectdirection='any',
#         legend=dict(
#             title='',
#             orientation='h',
#             yanchor="bottom",
#             y=1.02,
#             xanchor="right",
#             x=1
#         ),
#         xaxis_title="Completeness",
#         yaxis_title="Contamination"
#     )

#     # Distribution plot
#     if selected_data and selected_data['points']:
#         selected_indices = [point['pointIndex'] for point in selected_data['points']]
#         df_for_hist = plot_df.iloc[selected_indices]
#     else:
#         df_for_hist = plot_df

#     df_for_hist[color_var] = df_for_hist[color_var].replace(-999, 'Not available')
    
#     if is_numeric(df_for_hist[color_var]):
#         dist_fig = px.histogram(
#             df_for_hist[df_for_hist[color_var] != 'Not available'],
#             x=color_var,
#             title=f'Distribution of {color_var} for Selected Points',
#             nbins=30
#         )
#     else:
#         value_counts = df_for_hist[color_var].value_counts()
#         dist_fig = px.bar(
#             x=value_counts.index,
#             y=value_counts.values,
#             title=f'Distribution of {color_var} for Selected Points'
#         )
#         dist_fig.update_xaxes(title_text=color_var)
#         dist_fig.update_yaxes(title_text='Number of genomes')

#     not_available_count = (df_for_hist[color_var] == 'Not available').sum()
#     if not_available_count > 0:
#         dist_fig.add_annotation(
#             x=0.5, y=1.05,
#             xref='paper', yref='paper',
#             text=f'Number of genomes without metadata: {not_available_count}',
#             showarrow=False,
#             yshift=10,
#             font=dict(color='red', size=12)
#         )

#     dist_fig.update_layout(plot_bgcolor='white', font={'family': 'Arial', 'size': 14})

#     return scatter_fig, dist_fig    
 
@app.callback(
    Output('sunburst-plot', 'figure'),
    [Input('metadata-dropdown', 'value'),
     Input('reset-sunburst-button', 'n_clicks')]
)
def update_sunburst(selected_metadata, n_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return create_sunburst(df,selected_metadata)
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == 'reset-sunburst-button':
        return create_sunburst(df,selected_metadata)
    
    return create_sunburst(df,selected_metadata)

@app.callback(
    Output('sunburst-level-display', 'children'),
    [Input('sunburst-plot', 'clickData'),
     Input('reset-sunburst-button', 'n_clicks')]
)
def update_level_display(clickData, n_clicks):
    ctx = dash.callback_context
    if not ctx.triggered:
        return 'You are seeing now: Entire Dataset'
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == 'reset-sunburst-button':
        return 'You are seeing now: Entire Dataset'
    
    if clickData:
        current_path = clickData['points'][0]['id'].split('/')
        # Remove empty strings and 'None' values
        current_path = [p for p in current_path if p and p != 'None']
        if current_path:
            return f'You are seeing now: {" > ".join(current_path)}'
    
    return 'You are seeing now: Entire Dataset'

@app.callback(
    [Output('filtered-scatter-plot', 'figure'),
     Output('filtered-distribution-plot', 'figure')],
    [Input('sunburst-plot', 'clickData'),
     Input('metadata-dropdown', 'value'),
     Input('reset-sunburst-button', 'n_clicks')],
    [State('sunburst-plot', 'figure')]
)
def update_filtered_plots(clickData, selected_metadata, n_clicks, current_figure):
    ctx = dash.callback_context
    if not ctx.triggered:
        return dash.no_update, dash.no_update
    
    trigger_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if trigger_id == 'reset-sunburst-button':
        filtered_df = df.copy()
    else:
        filtered_df = df.copy()
        if clickData:
            current_path = clickData['points'][0]['id'].split('/')
            for i, level in enumerate(current_path):
                if level and i < len(tax_ranks):
                    filtered_df = filtered_df[filtered_df[tax_ranks[i]] == level]
    
    # if filtered_df.empty:
    #     return dash.no_update, dash.no_update
    
    if filtered_df.empty or selected_metadata is None:
        # Create empty figures with messages if there's no data or metadata
        scatter_fig = go.Figure()
        scatter_fig.add_annotation(text="No data to display or no metadata selected", 
                                  xref="paper", yref="paper",
                                  x=0.5, y=0.5, showarrow=False)
        
        dist_fig = go.Figure()
        dist_fig.add_annotation(text="No data to display or no metadata selected", 
                               xref="paper", yref="paper",
                               x=0.5, y=0.5, showarrow=False)
        
        return scatter_fig, dist_fig

    hover_data = ['formatted_classification'] + metadata_columns
    
    # Data preprocessing
    plot_df = filtered_df.copy()
    plot_df[selected_metadata] = plot_df[selected_metadata].replace('None', np.nan)
    is_numeric_metadata = pd.api.types.is_numeric_dtype(plot_df[selected_metadata].dropna())
    
    if is_numeric_metadata:
        # Numeric data processing
        plot_df[selected_metadata] = pd.to_numeric(plot_df[selected_metadata], errors='coerce')
        placeholder_mask = plot_df[selected_metadata].isna()
        
        # Main scatter plot
        scatter_fig = px.scatter(
            plot_df[~placeholder_mask],
            x='Completeness',
            y='Contamination',
            color=selected_metadata,
            hover_data=hover_data,
            title='Filtered: Completeness vs Contamination'
        )
        
        # Add a separate trace for NA values
        if placeholder_mask.any():
            scatter_fig.add_trace(
                go.Scatter(
                    x=plot_df[placeholder_mask]['Completeness'],
                    y=plot_df[placeholder_mask]['Contamination'],
                    mode='markers',
                    marker=dict(
                        color='red',
                        size=10,
                        symbol='x'
                    ),
                    name='Not available'
                )
            )
        
        # Distribution plot
        dist_fig = px.histogram(
            plot_df[~placeholder_mask],
            x=selected_metadata,
            title=f'Filtered Distribution of {selected_metadata}',
            nbins=30
        )
        
        na_count = placeholder_mask.sum()
        if na_count > 0:
            dist_fig.add_annotation(
                x=0.5,
                y=1.05,
                xref='paper',
                yref='paper',
                text=f'Number of genomes without metadata: {na_count}',
                showarrow=False,
                font=dict(color='red', size=12)
            )
    else:
        # String/categorical data processing
        scatter_fig = px.scatter(
            plot_df,
            x='Completeness',
            y='Contamination',
            color=selected_metadata,
            hover_data=hover_data,
            title='Filtered: Completeness vs Contamination'
        )
        
        # Categorical distribution plot
        dist_fig = px.histogram(
            plot_df,
            x=selected_metadata,
            title=f'Filtered Distribution of {selected_metadata}'
        )
    
    # Common layout setting
    scatter_fig.update_layout(
        plot_bgcolor='white',
        font={'family': 'Arial', 'size': 14},
        legend=dict(
            orientation='h',
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    dist_fig.update_layout(
        plot_bgcolor='white',
        font={'family': 'Arial', 'size': 14}
    )
    
    return scatter_fig, dist_fig


@app.callback(
    Output('data-table', 'data'),
    [Input('search-input', 'value'),
     Input('sunburst-plot', 'clickData'),
     Input('metadata-dropdown', 'value'),
     Input('data-table', 'sort_by')],
    [State('data-table', 'data')]
)
def update_table(search_value, clickData, selected_metadata, sort_by, current_data):
    filtered_df = df.copy()
    
    # Filtering based on Sunburst selection
    if clickData:
        current_path = clickData['points'][0]['id'].split('/')
        for i, level in enumerate(current_path):
            if level and i < len(tax_ranks):
                filtered_df = filtered_df[filtered_df[tax_ranks[i]] == level]
    
    # Filtering based on search value
    if search_value:
        filtered_df = filtered_df[filtered_df.apply(lambda row: any(str(search_value).lower() in str(cell).lower() for cell in row), axis=1)]
    
    # Apply sorting
    if sort_by:
        filtered_df = filtered_df.sort_values(
            [col['column_id'] for col in sort_by],
            ascending=[col['direction'] == 'asc' for col in sort_by],
            inplace=False
        )
    
    return filtered_df.to_dict('records')



@app.callback(
    Output('data-table', 'columns'),
    Input('column-selector', 'value')
)
def update_table_columns(selected_columns):
    valid_columns = [col for col in selected_columns if col != 'formatted_classification']
    return [{"name": col, "id": col} for col in valid_columns]
    #return [{"name": i, "id": i} for i in selected_columns if i != 'formatted_classification']

# download selected genome in the users local server callback
@app.callback(
    Output("download-dataframe-csv", "data"),
    Input("save-local-button", "n_clicks"),
    State("filename-input", "value"),
    State("data-table", "data"),
    prevent_initial_call=True,
)
def save_to_local(n_clicks, filename, table_data):
    if n_clicks == 0:
        raise PreventUpdate
    
    if not filename:
        filename = "data"

    filtered_df = pd.DataFrame(table_data)

    result_df = df_raw[df_raw['Genome'].isin(filtered_df['Genome'])]

    return dcc.send_data_frame(result_df.to_csv, f"{filename}.csv", index=False)

#    full_filtered_df = df.loc[filtered_df.index]

#    return dcc.send_data_frame(full_filtered_df.to_csv, f"{filename}.csv", index=False)

#    return dcc.send_data_frame(df.to_csv, f"{filename}.csv", index=False)

# server save callback original

# @app.callback(
#     Output("save-status", "children"),
#     Input("save-server-button", "n_clicks"),
#     State("filename-input", "value"),
#     State("data-table", "data"),
#     prevent_initial_call=True,
# )
# def save_to_server(n_clicks, filename, table_data):
#     if n_clicks == 0:
#         raise PreventUpdate
    
#     if not filename:
#         filename = "data"
    
#     df = pd.DataFrame(table_data)
#     server_path = os.path.join(os.getcwd(), f"{filename}.csv")
#     df = df.drop(columns=['formatted_classification'])
#     df.to_csv(server_path, index=False)
#     return f"File saved to server at {server_path}"

@app.callback(
    Output("overwrite-modal", "is_open"),
    [Input("save-server-button", "n_clicks"),
     Input("confirm-overwrite", "n_clicks"),
     Input("cancel-overwrite", "n_clicks")],
    [State("filename-input", "value"),
     State("overwrite-modal", "is_open")],
    prevent_initial_call=True
)
def toggle_modal(save_clicks, confirm_clicks, cancel_clicks, filename, is_open):
    ctx = dash.callback_context
    if not ctx.triggered:
        return is_open
        
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if button_id == "save-server-button":
        if not filename:
            filename = "genome_selector_result"
        
        server_path = os.path.join(os.getcwd(), f"{filename}.csv")
        return os.path.exists(server_path)
    
    elif button_id == "confirm-overwrite" or button_id == "cancel-overwrite":
        return False
    
    return is_open


@app.callback(
    Output("save-status", "children"),
    [Input("save-server-button", "n_clicks"),
     Input("confirm-overwrite", "n_clicks"),
     Input("cancel-overwrite", "n_clicks")],
    [State("filename-input", "value"),
     State("data-table", "data"),
     State("overwrite-modal", "is_open")],
    prevent_initial_call=True
)
def save_to_server(save_clicks, confirm_clicks, cancel_clicks, filename, table_data, modal_is_open):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate
        
    button_id = ctx.triggered[0]['prop_id'].split('.')[0]
    
    if not filename:
        filename = "genome_selector_result"
    
    server_path = os.path.join(os.getcwd(), f"{filename}.csv")
    
    if button_id == "save-server-button" and not os.path.exists(server_path):
        # If the file does not exist, save it immediately
        filtered_df = pd.DataFrame(table_data)
        
        # Remove the formatted_classification column
        if 'formatted_classification' in filtered_df.columns:
            filtered_df = filtered_df.drop(columns=['formatted_classification'])
        
        # Extract rows matching the original data
        result_df = df_raw[df_raw['Genome'].isin(filtered_df['Genome'])].copy()
        
        # Remove the formatted_classification column from the original data
        if 'formatted_classification' in result_df.columns:
            result_df = result_df.drop(columns=['formatted_classification'])
        
        result_df.to_csv(server_path, index=False)
        return f"File saved to server at {server_path}"
    
    elif button_id == "confirm-overwrite" and not modal_is_open:
        # Confirm overwrite and save
        filtered_df = pd.DataFrame(table_data)
        
        # Remove the formatted_classification column
        if 'formatted_classification' in filtered_df.columns:
            filtered_df = filtered_df.drop(columns=['formatted_classification'])
        
        # Extract rows matching the original data
        result_df = df_raw[df_raw['Genome'].isin(filtered_df['Genome'])].copy()
        
        # Remove the formatted_classification column from the original data
        if 'formatted_classification' in result_df.columns:
            result_df = result_df.drop(columns=['formatted_classification'])
        
        result_df.to_csv(server_path, index=False)
        return f"File overwritten at {server_path}"
    
    elif button_id == "cancel-overwrite":
        return "Save operation cancelled"
    
    return dash.no_update


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='metaFun: genome selector for COMPARATIVE_ANNOTATION')
    parser = argparse.ArgumentParser(description='refer to the documentation at https://metafun-doc-v01.readthedocs.io/en/latest/')
    parser.add_argument('-i', '--input', help='Input CSV file', required=True)
    parser.add_argument('-p', '--port', help='Port to run the server on', type=int, default=8050)

    args = parser.parse_args()
    app.run_server(debug=True, host='0.0.0.0', port=args.port)
    