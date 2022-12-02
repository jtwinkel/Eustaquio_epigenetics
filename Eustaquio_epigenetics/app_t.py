#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 15:58:01 2022

@author: jonwinkelman
"""



import pysam
from plotly import graph_objects as go
import plotly.offline as pyo
import pandas as pd
from Bio import SeqIO
import dash
from dash import dcc, html
import plotly.graph_objs as go
from dash.dependencies import Input, Output, State


path_to_genbank = './data/references/Reference_FERM_BP3421.gbk'
bam_filepath = './RNAseq_bam_files/BAN-1.bam'

def pileup_on_each_contig(bam_filepath):
    """
    return a df with number of reads at each position in each chromosome
    df structure: 
        index: genomic postion (one-based numbering)
        columns: each column represents one chromosome, plasmid, etc
        
        index starts at 0, the first nucleotide, and ends at the last nucleotide
        of the longest chromosome, plasmid, etc
        
        if there are no reads a genomic position within the length of the contig,
        then the reads will be equal to 0. If the genomic postion is out of the
        range of contig, then read value at that position will be equal to np.nan
    

    """
    ban1_file = pysam.AlignmentFile(bam_filepath, 'rb')
    lst_of_dfs = []
    
    for contig in ban1_file.references:
        pileup_genome = {}
        for col in ban1_file.pileup(contig = contig):
            if pileup_genome.get(contig):
                pileup_genome[contig].update({col.reference_pos:col.nsegments})
            else:
                pileup_genome[contig] = {col.reference_pos:col.nsegments}
        chromosome_length = ban1_file.get_reference_length(contig)
        df = pd.DataFrame(pileup_genome)
        df = df.reindex(index=range(chromosome_length), fill_value=0)
        lst_of_dfs.append(df)
        
    df = pd.concat(lst_of_dfs, axis=1, join='outer')
    df.index = df.index +1 #change to 1-based numbering
    df.index.names = ['position']
    df.to_csv('./data/pilup_each_chromosome.csv')
    return df





# =============================================================================
# Non-callback functions
# ============================================================================= 
       
# =============================================================================
# df = pd.DataFrame.from_dict(pileup_genome, orient='index')
# start = 0
# stop = 10000
# data = go.Bar(
#         x = df.index[start:stop],
#         y = df.iloc[start:stop,0]
#     )
# 
# fig = go.Figure(data = data)
# #pyo.plot(fig)
# =============================================================================

class make_chromosome_dict():
    
    def __init__(self, path_to_genbank):    
        self.path_to_genbank = path_to_genbank
        self.full_dict = self.build_dict()
        self.chromosomes = list(self.full_dict.keys())
        
    def build_dict(self):
        """
        builds dict {chromosomeID:{featureID:feature}} from a genbank annotation file
        
        """
        chromosome_dict = {}
        for seq_record in SeqIO.parse(self.path_to_genbank , "genbank"):  
            for feature in seq_record.features:
                if feature.type == 'gene':
                    if chromosome_dict.get(seq_record.name): 
                        chromosome_dict[seq_record.name].update({feature.qualifiers['locus_tag'][0]:feature})
                    else:
                        chromosome_dict[seq_record.name] = {feature.qualifiers['locus_tag'][0]:feature}
        return chromosome_dict


    
    



# =============================================================================
# Main app
# =============================================================================           
    
   
app = dash.Dash(__name__)
assembly_obj = make_chromosome_dict(path_to_genbank)
chromosome_dict = assembly_obj.full_dict


chromosome_gene_selector = html.Div([
                           html.H4('Chromosome selector'),
                           dcc.Dropdown(
                                id='Chromosome_selector',
                                options=[{'label':chromosome, 'value':chromosome} for chromosome in chromosome_dict.keys()],
                                value = list(chromosome_dict.keys())[0]
                                ),
                                ])
                        
feature_selector = html.Div([
                        html.H4('Feature list'),      
                        dcc.Dropdown(
                            id='opt_feature_dropdown',
                            ),
                            ]
                        )

chromosome_map = html.Div([
                    html.H4('Chromosome pilup'),
                    dcc.Graph(id = 'chromosome_pileup')
                    ]
                    )

                                
app.layout  = html.Div([
    chromosome_gene_selector,
    feature_selector,
    chromosome_map
    
    
    ]) #, style={
         #    'display': 'inline-block',
          #   'padding-left': 50})
                                     

# =============================================================================
# callback functions
# =============================================================================                        
                    
@app.callback(
    Output('opt_feature_dropdown', 'options'),
    [Input('Chromosome_selector', 'value')]
)
def set_features_for_dropdown(chromosome_ID):
    'return list of features from selected chromosome for feature dropdown'
    return [{'label': feature.strip(), 'value': feature.strip()} for feature in chromosome_dict[chromosome_ID].keys()]



@app.callback(
    Output('chromosome_pileup', 'figure'),
    [Input('Chromosome_selector', 'value'),
     Input('opt_feature_dropdown', 'value')
     ]
)
def get_feature_for_map(chromosome_ID, feature_ID):
    'return list of features from selected chromosome for feature dropdown'
    print(chromosome_ID, feature_ID)
    offset = 200
    df = pd.read_csv('./data/pilup_each_chromosome.csv').set_index('position')
    feature = chromosome_dict[chromosome_ID][feature_ID]
    start = feature.location.start - offset
    end = feature.location.end + offset
    strand = feature.location.strand
    data = go.Bar(
             x = df.index[start:end],
             y = df.loc[start:end,chromosome_ID]
         )
    layout = go.Layout(title = f'chromosome: {chromosome_ID} feature ID: {feature_ID}')
    fig = go.Figure(data = data, layout = layout)
    return fig


if __name__ == '__main__':
    app.run_server()