#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  3 08:38:52 2022

@author: jonwinkelman
"""
import bisect
import dash
import json
from dash import dcc, html, callback_context as ctx
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State
import pandas as pd
import numpy as np
import os
import zipfile
import tempfile
from jw_utils import dna_utils as du
from jw_utils import parse_fasta as pf
from jw_utils import file_utils as fu
from jw_utils import parse_gff as pgf
from jw_utils import ribo_profile3 as rp
from jw_utils import plotly_utils as plu
from jw_utils import genome_utils as gu

from jw_utils import app_functions as afns
from jw_utils import app_elements as ae

# =============================================================================
# Add the correct path to "path_to_results"
# set lab_name variable
# =============================================================================
#path_to_results = './data/Orthofinder_data/Proteomes/OrthoFinder/Results_Apr28_1'
lab_name = 'Mankin Lab'
USERNAME_PASSWORD_PAIRS = [['Trestle', 'MankinLab']]
# =============================================================================
# |^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|^|
# =============================================================================

path_to_references = './data/references'
path_to_valcounts_dir = './data/value_counts'
path_to_termination_dir = './data/termination_files'
path_to_gff = os.path.join(path_to_references, 'GCF_013166975.1_ASM1316697v1_genomic.gff')
path_to_proteome = os.path.join(path_to_references, 'GCF_013166975.1.faa')
path_to_fasta_genome = os.path.join(path_to_references, 'BL21-CP053601.1.fa')
list_of_valcount_paths = fu.get_filepaths_in_dir(path_to_valcounts_dir, end_to_exclude='.DS_Store')

colors = {
't_blue': 'rgba(0,102,153,255)',
't_green': 'rgba(61,174,43,255)',
'orange':'#F9B257',
't_red': 'rgb(255,20,20)',
'seagreen':'#2c8d42',

}




# =============================================================================
# functions
# =============================================================================



# =============================================================================
# function calls - be careful, these are global variables???c
# =============================================================================

rp_object_dict = afns.rb_obj_store(list_of_valcount_paths, path_to_gff).rp_object_dict
afns.check_for_proper_files(rp_object_dict,path_to_termination_dir, 
                            path_to_proteome, path_to_fasta_genome)
    
    
density_filter_text = 'number of reads per nucleotide, per million total reads'
neighbor_filter_text = ('the minimum intergene distance allowed.\n'
                        'note, this distance is negative if genes overlap\n' 
                        'on the same strand')
# =============================================================================
#                               Main app
# =============================================================================

app = dash.Dash(__name__)
#auth =dash_auth.BasicAuth(app, USERNAME_PASSWORD_PAIRS)
server = app.server

header_container = ae.create_header('Mankin Lab')
feature_dropdown = ae.create_feature_dropdown(path_to_gff)
up_downstream_input = ae.create_up_downstream_input()
experiment_selector = ae.create_experiment_selector(path_to_valcounts_dir)
# =============================================================================
#                               Main app layout
# =============================================================================

app.layout = html.Div([
    header_container,
    html.Div([      
        html.Button("Download data", id="btn-download-txt"), #styled in external css file
        dcc.Download(id="download_termscores1"),
        dcc.Download(id="download_termscores2"),
        dcc.Download(id="download_start_aligned1"),
        dcc.Download(id="download_start_aligned2"),
        dcc.Download(id="download_start_aligned3"),
        dcc.Download(id="download_start_aligned4"),
        dcc.Download(id="download_stop_aligned1"),
        dcc.Download(id="download_stop_aligned2"),
        dcc.Download(id="download_stop_aligned3"),
        dcc.Download(id="download_stop_aligned4"),
        ], style = {'width': '100%',
                    'display':'block',
                    'margin-left':'auto',
                    'margin-right':'auto'}
    ),
    html.H4('Term scores density_filter, (RPMs)', style={'textAlign': 'left'}),
    html.H5(density_filter_text, style={'textAlign': 'left'}),
    dcc.Input(id='ts_graph_density_filter',
        type='number',
        value = 0),
    html.H4('Term scores neighbor filter, (nts)', style={'textAlign': 'left'}),
    html.H5(neighbor_filter_text, style={'textAlign': 'left'}),
    dcc.Input(id='ts_neighbor_filter',
        type='number',
        value = -5000), 
    dcc.Graph(id='term_scores', figure=afns.make_term_score_fig(path_to_termination_dir, path_to_gff)),
    feature_dropdown,
    up_downstream_input,
    experiment_selector,
    dcc.Graph(id='gene_profile_graph'),
    html.H4('Density_filter, (RPMs)', style={'textAlign': 'left'}),
    dcc.Input(id='graph_density_filter',
        type='number',
        value = 0),
    html.H4('Neighbor_filter, (nucleotides)', style={'textAlign': 'left'}),
    dcc.Input(id='neighbor_filter',
        type='number',
        value = -5000), 
    html.Div([
        html.Div([
            dcc.Graph(id='start_graph'),
        ], style = {'width':'70%',
                    'display':'inline-block',}),
        html.Div([
            dcc.Graph(id='start_position_dist')
        ], style = {'width':'15%',
                    'display':'inline-block',}),
    ],style = {'width':'100%'}),        
    html.Div([
        html.Div([
            dcc.Graph(id='stop_graph')
        ], style = {
            'width':'70%',
            'display':'inline-block',}),
        html.Div([
            dcc.Graph(id='stop_position_dist')
        ], style = {
            'width':'15%',
            'display':'inline-block',}),
    ],style = {'width':'100%'}),
    html.Div(id='test'),
    html.Div(id='test_div'),
    dcc.Graph(id='gene_profile_graph2'),
], style = {'padding':25})


# =============================================================================
# callback functions
# =============================================================================

@app.callback(
    Output('term_scores', 'figure'),
    [Input('ts_graph_density_filter','value'),
     Input('ts_neighbor_filter','value')
     ])
def update_termscore_graph(filt_value, clearance):
    df_dict, fig =afns.make_term_score_fig(path_to_termination_dir,
                                 path_to_gff,
                                 filt_value=filt_value,
                                 clearance=clearance,
                                 return_type = 'all')
    arrs = ('x_min', 'x_max', 'y_min', 'y_max')
    d = {ar:[] for ar in arrs}
    for name, df in df_dict.items():
        d['x_min'].append(  np.log2(df.iloc[:,0].min() )  )
        d['x_max'].append(  np.log2(df.iloc[:,0].max() ) )
        d['y_min'].append(  np.log2(df.iloc[:,2].min() ) )
        d['y_max'].append(  np.log2(df.iloc[:,2].max() ) )
    fig.add_trace(go.Scatter(
        x = [min(d['x_min']), max(d['x_max'])],
        y = [min(d['y_min']), max(d['y_max'])],
        mode = 'lines',
        line = dict(color='rgba(100,100,100,0.1)'),
        name = 'Diagonal'
        ))
    return fig



@app.callback(
    Output('stop_graph', 'figure'),
    [Input('graph_density_filter','value'),
     Input('neighbor_filter','value')
     ])
def update_stop_graph(filt_value, clearance):
    filt_value = float(filt_value)
    if filt_value <=0: density_filt=False
    else: density_filt=True
    clearance = int(clearance)
    fig =  afns.plot_stops_starts_metagene(path_to_gff,  stops_or_start='stops',
                                    density_filt=density_filt,
                                    filt_level=filt_value,
                                    neighbor_clearance=clearance
                                    )
    return afns.update_fig_style(fig, 
                            title='Stop-aligned genes', 
                            ylabel='Normalized reads',
                            xlabel='Position relative to stop')




@app.callback(
    Output('start_graph', 'figure'),
    [Input('graph_density_filter','value'),
     Input('neighbor_filter','value')
     ])
def update_start_graph(filt_value, clearance):
    filt_value = float(filt_value)
    if filt_value <=0: density_filt=False 
    else: density_filt=True    
    clearance = int(clearance)
    fig =  afns.plot_stops_starts_metagene(path_to_gff,  stops_or_start='starts',
                               density_filt=density_filt,
                               filt_level=filt_value,
                               neighbor_clearance=clearance
                               )
    return afns.update_fig_style(fig, 
                            title='Start-aligned genes', 
                            ylabel='Normalized reads',
                            xlabel='Position relative to start')




@app.callback(
    Output('gene_profile_graph2', 'figure'),
    [Input('start_position_dist','clickData'),
     Input('stop_position_dist','clickData'),
     Input('up_input','value'),
     Input('down_input','value')
     ],
    )
def graph_gene_region2(start_clickData, stop_clickData, up_input, down_input):
    '''
    update the second feature graph with clickdata from the start or stop 
    position distribution graph
    '''
    gff_obj=pgf.make_seq_object_dict(path_to_gff, feature_type='CDS')
    clickData = None
    fig = go.Figure()
    if ctx.triggered[0]['prop_id'] == 'start_position_dist.clickData':
        clickData = start_clickData
        graph_type = 'start'
    if ctx.triggered[0]['prop_id'] == 'stop_position_dist.clickData':
        clickData = stop_clickData
        graph_type = 'stop'
        
    if clickData:
        #print(clickData)
        t_lst = clickData['points'][0]['text'].split('<br>') 
        feature_id = t_lst[0]
        experiment = t_lst[1]
        selected_position = int(t_lst[2].split(' ')[-1])
        start_offset, stop_offset, trimmed_df = afns.get_offsets(feature_id, up_input, down_input, path_to_gff)
        l = []
        df, start, stop = rp_object_dict[experiment].gene_hits(feature_id,
                                 reads_per_mil=True,start_offset=start_offset,
                                 stop_offset=stop_offset)
        fig.add_trace(go.Bar(
            x=df.index,
            y=df.loc[:,df.columns[0]],
            name = f'{experiment}\n{feature_id}'
            ))
        l.append(df.loc[:,df.columns[0]].max())
        max_val = max(l)
        if trimmed_df.loc[feature_id, 'strand'] == '+':
            if graph_type == 'start':
                position_line =trimmed_df.loc[feature_id, 'start'] + selected_position
            elif graph_type == 'stop':
                position_line =trimmed_df.loc[feature_id, 'end'] + selected_position
        elif trimmed_df.loc[feature_id, 'strand'] == '-':
            if graph_type == 'start':
                position_line =trimmed_df.loc[feature_id, 'end'] - selected_position
            elif graph_type == 'stop':
                position_line =trimmed_df.loc[feature_id, 'start'] - selected_position
                
        fig.add_trace(go.Scatter(
            x=[position_line,position_line],
            y=[0,max_val*1.2],
            name = f'Selected position: {selected_position}, {graph_type}-aligned',
            mode = 'lines',
            line = {'width': 5,
                    'color':'rgba(100,100,100,0.5)'}
            ))
                 
        fig.add_trace(go.Scatter(
            x=[stop,stop],
            y=[0,max_val*1.2],
            name = 'stop',
            mode = 'lines',
            line = {'width': 3,
                    'dash':'dash',
                    'color':'rgba(200,100,100,0.25)'}
            ))
        fig.add_trace(go.Scatter(
            x=[start,start],
            y=[0,max_val*1.2],
            name = 'start',
            mode = 'lines',
            line = {'width': 3,
                    'dash':'dash',
                    'color':'rgba(100,200,100,0.5)'}
            ))
        #plot the feature arrows
        for ind in trimmed_df.index:
            if trimmed_df.loc[ind,'strand'] == '-':
                stop = trimmed_df.loc[ind,'start']
                start = trimmed_df.loc[ind,'end']
            else:
                start = trimmed_df.loc[ind,'start']
                stop = trimmed_df.loc[ind,'end']   
            fig.add_trace(go.Scatter(
                x = plu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[0],
                y = plu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[1],
                fill="toself",
                mode = 'lines',
                line = {'width': 1.5,
                        'color':'rgb(200,200,200)',
                        },
                name=ind,
                showlegend=False,
                text = f'{gff_obj[ind].ID}<br>{gff_obj[ind].locus_tag}<br>{gff_obj[ind].Name}<br>{gff_obj[ind].product}'
                ))
        return afns.update_fig_style(fig, 
                        title=ctx.triggered[0]['prop_id'], 
                        ylabel='Reads (per million',
                        xlabel='Genomic position')
    else:
        return afns.update_fig_style(fig, 
                        title='Selected gene region',
                        ylabel='Reads (per million',
                        xlabel='Genomic position')





@app.callback(
        Output('start_position_dist', 'figure'),
        [Input('start_graph','clickData'),
         Input('graph_density_filter','value'),
         Input('neighbor_filter','value')]
        )
def update_start_distribution(clickData, filt_level, clearance):
    '''
    update distribution of reads for each gene, generated from a selected
    genomic position on the stop-aligned graph
    '''
    starts = './data/metagene_stops_starts/starts'
    exp_dict = afns.simple_experimentalSetup_dict(starts)
    filt_level = float(filt_level)
    filt_gene_dict = afns.density_filter(filt_level)
    genes_clearance = gu.get_genes_wo_neighbors(path_to_gff, clearance=clearance)
    fig = go.Figure()
    if clickData:
        position = clickData['points'][0]['x']
        exp = clickData['points'][0]['curveNumber']
        names = []
        for condition in exp_dict.keys():  
            for replicate in exp_dict[condition]:
                names.append(afns.expname_from_path(replicate)) 
        selected_exper = {pos:ex for pos,ex in zip(range(len(names)),names)}[exp]
        genes_to_keep = set(filt_gene_dict[selected_exper]).intersection(set(genes_clearance))
        path = fu.get_filepaths_in_dir(starts, begin_to_keep=selected_exper)[0]
        df = pd.read_csv(path)
        ser=df.set_index(df.columns[0]).loc[:,genes_to_keep].loc[position,:]
        text = [f'{feature}<br>{selected_exper}<br>selected position: {position}' for feature in ser.keys()]
        fig.add_trace(go.Scatter(
            x = np.random.normal(loc=position, scale=0.2, size = ser.shape[0]),
            y = ser,
            text = text,
            mode = 'markers',
            marker=dict(color='rgba(100,100,100,0.3)')
            )
        )
        fig = afns.update_fig_style(fig, 
                            title=f'Reads at selected<br>position {position}', 
                            ylabel='Normalized reads',
                            xlabel='Pos. relative to start')
        return fig
    else:
        fig.add_trace(go.Scatter(
            x = [],
            y = [],
            mode='markers')
            )
        return afns.update_fig_style(fig, 
                            title='Reads at selected<br>position', 
                            ylabel='Normalized reads',
                            xlabel='Pos. relative to start')      
        
        
@app.callback(
        Output('stop_position_dist', 'figure'),
        [Input('stop_graph','clickData'),
         Input('graph_density_filter','value'),
         Input('neighbor_filter','value')
         ]
        )
def update_stop_distribution(clickData, filt_level, clearance):
    '''
    update distribution of reads for each gene, generated from a selected
    genomic position on the stop-aligned graph
    '''
    stops = './data/metagene_stops_starts/stops'
    exp_dict = afns.simple_experimentalSetup_dict(stops)
    filt_level = float(filt_level)
    filt_gene_dict = afns.density_filter(filt_level)
    genes_clearance = gu.get_genes_wo_neighbors(path_to_gff, clearance=clearance)
    if clickData:
        fig = go.Figure()
        position = clickData['points'][0]['x']
        exp = clickData['points'][0]['curveNumber']
        names = []
        for condition in exp_dict.keys():  
            for replicate in exp_dict[condition]:
                names.append(afns.expname_from_path(replicate))
                
        selected_exper = {pos:ex for pos,ex in zip(range(len(names)),names)}[exp]
        genes_to_keep = set(filt_gene_dict[selected_exper]).intersection(set(genes_clearance))
        path = fu.get_filepaths_in_dir(stops, begin_to_keep=selected_exper)[0]
        df = pd.read_csv(path)
        ser=df.set_index(df.columns[0]).loc[:,genes_to_keep].loc[position,:]
        text = [f'{feature}<br>{selected_exper}<br>selected position: {position}' for feature in ser.keys()]
        fig.add_trace(go.Scatter(
            x = np.random.normal(loc=position, scale=0.2, size = ser.shape[0]),
            y = ser,
            text = text,
            mode = 'markers',
            marker=dict(color='rgba(100,100,100,0.3)')
            )
        )
        return afns.update_fig_style(fig, 
                            title=f'Reads at selected<br>position {position}', 
                            ylabel='Normalized reads',
                            xlabel='Pos. relative to stop')
         
    else:
        fig = go.Figure()
        fig.add_trace(go.Scatter(
            x = [],
            y = [],
            mode='markers')
            )
        return afns.update_fig_style(fig, 
                            title='Reads at selected<br>position', 
                            ylabel='Normalized reads',
                            xlabel='Pos. relative to stop')

    
@app.callback(
        Output('feature_dropdown', 'value'),
        [Input('term_scores','clickData')]
        )
def update_dropdown_w_termscore_feature(clickData):
    'clickData {points:[{},]'
    if clickData: 
        t_lst = clickData['points'][0]['text'].split('<br>')
        value = t_lst[1]
        return value


@app.callback(
        Output('experiment_selector', 'value'),
        #Output('test_div', 'children'),
        [Input('term_scores','clickData')]
        )
def update_dropdown_w_termscore_exp(clickData):
    'clickData {points:[{},]'
    if clickData: 
        t_lst = clickData['points'][0]['text'].split('<br>')
        t = t_lst[0].replace('_','-')
        value = [f'{t}-1', f'{t}-2']
        return value


@app.callback(
    Output('gene_profile_graph', 'figure'),
    [Input('experiment_selector','value'),
     Input('feature_dropdown','value'),
     Input('up_input','value'),
     Input('down_input','value')
     ],
    )
def graph_gene_region(experiment_name, feature_id, up_input, down_input):
    fig = go.Figure()
    #if up_input or down_input:     
    gff_obj=pgf.make_seq_object_dict(path_to_gff, feature_type='CDS')
    if experiment_name and feature_id:
        start_offset, stop_offset, trimmed_df = afns.get_offsets(feature_id, up_input, down_input, path_to_gff)
        l = []
        for name in experiment_name:
            
            df, start, stop = rp_object_dict[name].gene_hits(feature_id,
                                     reads_per_mil=True,start_offset=start_offset,
                                     stop_offset=stop_offset)
            fig.add_trace(go.Bar(
                x=df.index,
                y=df.loc[:,df.columns[0]],
                name = f'{name}\n{feature_id}'
                ))
            l.append(df.loc[:,df.columns[0]].max())
            max_val = max(l)
        fig.add_trace(go.Scatter(
            x=[stop,stop],
            y=[0,max_val*1.2],
            name = 'stop',
            mode = 'lines',
            line = {'width': 3,
                    'dash':'dash',
                    'color':'rgba(200,100,100,0.5)'}
            ))
        fig.add_trace(go.Scatter(
            x=[start,start],
            y=[0,max_val*1.2],
            name = 'start',
            mode = 'lines',
            line = {'width': 3,
                    'dash':'dash',
                    'color':'rgba(100,200,100,0.5)'}
            ))
        for ind in trimmed_df.index:
            if trimmed_df.loc[ind,'strand'] == '-':
                stop = trimmed_df.loc[ind,'start']
                start = trimmed_df.loc[ind,'end']
            else:
                start = trimmed_df.loc[ind,'start']
                stop = trimmed_df.loc[ind,'end']   
            fig.add_trace(go.Scatter(
                x = plu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[0],
                y = plu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[1],
                fill="toself",
                mode = 'lines',
                line = {'width': 1.5,
                        'color':'rgb(200,200,200)',
                        },
                name=feature_id,
                showlegend=False,
                text = f'{gff_obj[ind].ID}<br>{gff_obj[ind].locus_tag}<br>{gff_obj[ind].Name}<br>{gff_obj[ind].product}'
                ))
        fig = afns.update_fig_style(fig,
                   title='Gene (selected from termination score graph)',
                   ylabel = "Genomic position",
                   xlabel = "Reads (per million)"
                   )
    return fig




# =============================================================================
# DOWNLOAD DATA  clumsy...callback for each csv to download
# =============================================================================

@app.callback(
    Output("download_termscores1", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_termscores1(n_clicks, filt_value, clearance):
    df_dict = afns.make_term_score_fig(path_to_termination_dir,
                                 path_to_gff,
                                 filt_value=filt_value,
                                 clearance=clearance,
                                 return_type = 'df')
    condition = 'S1-Mono'
    df = df_dict[condition]
    return dcc.send_data_frame(df.to_csv, f"{condition}_termscores_densfilt_{filt_value}_clear_{clearance}.csv")


@app.callback(
    Output("download_termscores2", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_termscores2(n_clicks, filt_value, clearance):
    df_dict = afns.make_term_score_fig(path_to_termination_dir,
                                 path_to_gff,
                                 filt_value=filt_value,
                                 clearance=clearance,
                                 return_type = 'df')
   
    condition = 'S2-Disome'
    df = df_dict[condition]
    return dcc.send_data_frame(df.to_csv, f"{condition}_termscores_densfilt_{filt_value}_clear_{clearance}.csv")



@app.callback(
    Output("download_start_aligned1", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_start_aligned1(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='starts',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')

    experiment = 'S1-Mono-1'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_start_densfilt_{filt_value}_clear_{clearance}.csv")


@app.callback(
    Output("download_start_aligned2", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_start_aligned2(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='starts',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')
    experiment='S1-Mono-2'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_start_densfilt_{filt_value}_clear_{clearance}.csv")



@app.callback(
    Output("download_start_aligned3", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_start_aligned3(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='starts',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')
    experiment='S2-Disome-1'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_start_densfilt_{filt_value}_clear_{clearance}.csv")

@app.callback(
    Output("download_start_aligned4", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_start_aligned4(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='stops',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')
    experiment='S2-Disome-2'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_start_densfilt_{filt_value}_clear_{clearance}.csv")

@app.callback(
    Output("download_stop_aligned1", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_stop_aligned1(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='stops',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')

    experiment = 'S1-Mono-1'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_stop_densfilt_{filt_value}_clear_{clearance}.csv")


@app.callback(
    Output("download_stop_aligned2", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_stop_aligned2(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='stops',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')
    experiment='S1-Mono-2'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_stop_densfilt_{filt_value}_clear_{clearance}.csv")



@app.callback(
    Output("download_stop_aligned3", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_stop_aligned3(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='stops',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')
    experiment='S2-Disome-1'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_stop_densfilt_{filt_value}_clear_{clearance}.csv")

@app.callback(
    Output("download_stop_aligned4", "data"),
    [Input("btn-download-txt", "n_clicks")],
    [State('graph_density_filter','value'),
     State('neighbor_filter','value')
     ],
    prevent_initial_call=True,
)
def download_stop_aligned4(n_clicks, filt_value, clearance):
    start_df_dict = afns.plot_stops_starts_metagene(path_to_gff,
                                         stops_or_start='stops',
                                         density_filt=True, 
                                         filt_level=filt_value,
                                         neighbor_clearance=clearance,
                                         return_type='df')
    experiment='S2-Disome-2'
    df = start_df_dict[experiment]
    return dcc.send_data_frame(df.to_csv, f"{experiment}_stop_densfilt_{filt_value}_clear_{clearance}.csv")



if __name__ == '__main__':
    app.run_server()
