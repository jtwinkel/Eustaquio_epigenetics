#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 10 19:12:13 2022

@author: jonwinkelman
"""
import pandas as pd
from jw_utils import gene_profile as gpro
from plotly import offline as pyo
from plotly import graph_objects as go
import numpy as np
import os
from jw_utils import app_functions as afns
from jw_utils import parse_gff as pgf
from jw_utils import parse_fasta as pf
from jw_utils import genome_utils as gu
from jw_utils import parse_gbk as pgb
from jw_utils import plotly_utils as pu
from jw_utils import file_utils as fu
from jw_utils import gene_profile as gp
from Bio import SeqIO
import bisect

contigs=['BF000000.1', 'BF000000.2', 'BF000000.3', 'BF000000.4', 'BF000000.5']
path_to_valcount_dir = './results/analysis/value_counts'
    
def merge_exper_valcounts(path_to_valcount_dir, contig_name):
    """Return a merged df with val_cnts from each experiment for the given contig."""
    vc_fpaths = fu.get_filepaths_in_dir(path_to_valcount_dir, end_to_keep='.csv')
    i=0
    for path in vc_fpaths:
        if path.split('/')[-1].startswith(contig_name):
            print(path)
            i+=1
            suf = '_'.join(path.split('/')[-1].split('_')[1:3])
            if i<=1:
                df = pd.read_csv(path)
                df = df.set_index(df.columns[0])
                df.columns = [col+'_'+suf for col in df.columns]
            else:
                t=pd.read_csv(path)
                t = t.set_index(t.columns[0])
                t.columns = [col+'_'+suf for col in t.columns]
                df = pd.merge(df, t, left_index=True, right_index=True,  how='outer')
    return df

def get_read_sums(path_to_valcount_dir, contig_name):
    """Return sum of reads for each genomic position on + and - strands from merged df.
    
    sums reads at each position across all pathnames in the dir that start with the contig
    name. I.e. it sums across all experiments on a given contig.
    """
    merged_df = merge_exper_valcounts(path_to_valcount_dir, contig_name)
    df_t = pd.DataFrame()
    df = pd.DataFrame()
    df_b = pd.DataFrame()
    for col in merged_df.columns:
        if col.startswith("5'_(+)"):
            df_t[col] = merged_df[col]
        if col.startswith("5'_(-)"):
            df_b[col] = merged_df[col]
    df["5'_(+)_sum"] =df_t.sum(axis=1)
    df["5'_(-)_sum"] = df_b.sum(axis=1)
    return df

contig_name = 'BF000000.1'
def graph_gene_region_sum(feature_id, up_input, down_input, path_annot_file,
                          path_to_valcount_dir, contig_name=None, 
                          annot_type='gbk', feature_type = 'CDS'):
    """Plot the reads in region of the input feature with upstream and downstream genes.
    
    Arguments:
        df (pd.DataFrame): index should be the genomic position, each column 
        should be the value counts that we want mapped
        
    """
    df = get_read_sums(path_to_valcount_dir, contig_name)
    #create annotation object for contig
    gb_obj = SeqIO.parse(path_annot_file, "genbank")
    gb_contig = [contig for contig in gb_obj if contig.name==contig_name][0]
    #get feature annots
    feature = [feat for feat in gb_contig.features if feat.qualifiers.get('locus_tag')==[feature_id] and feat.type==feature_type][0]
    start = int(feature.location.start)
    end = int(feature.location.end)
    if start>end:
        start_codon=end
        stop_codon=start
    else:
        start_codon=start
        stop_codon=end
    
    #For generating annotation arrows etc
    start_offset, stop_offset, trimmed_df = afns.get_offsets(feature_id, up_input, 
                                        down_input, path_annot_file, annot_type=annot_type,
                                        contig_name=contig_name)
    l = []
    fig = go.Figure()
    #Actual read data
    df = get_read_sums(path_to_valcount_dir, contig_name)
    start_index = bisect.bisect_left(df.index,(start_codon - start_offset))
    end_index = bisect.bisect_right(df.index, (stop_codon + stop_offset))
    df = df.iloc[start_index:end_index,:]
    fig.add_trace(go.Bar(
        x=df.index,
        y=df.loc[:,df.columns[0]],
        name = f'{contig_name}\n{feature_id}'
        ))
    l.append(df.loc[:,df.columns[0]].max())  
    max_val = max(l)
    #plot stop and start codon lines on graph
    fig.add_trace(go.Scatter(
        x=[stop_codon,stop_codon],
        y=[0,max_val*1.2],
        name = 'stop',
        mode = 'lines',
        line = {'width': 3,
                'dash':'dash',
                'color':'rgba(200,100,100,0.5)'}
        ))
    fig.add_trace(go.Scatter(
        x=[start_codon,start_codon],
        y=[0,max_val*1.2],
        name = 'start',
        mode = 'lines',
        line = {'width': 3,
                'dash':'dash',
                'color':'rgba(100,200,100,0.5)'}
        ))
    #plot annoation arrows on graph using trimmed df
    for ind in trimmed_df.index:
        if trimmed_df.loc[ind,'strand'] == '-':
            stop = trimmed_df.loc[ind,'start']
            start = trimmed_df.loc[ind,'end']
        else:
            start = trimmed_df.loc[ind,'start']
            stop = trimmed_df.loc[ind,'end']   
        fig.add_trace(go.Scatter(
            x = pu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[0],
            y = pu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[1],
            fill="toself",
            mode = 'lines',
            line = {'width': 1.5,
                    'color':'rgb(200,200,200)',
                    },
            name=feature_id,
            showlegend=False,
            #text = f'{gff_obj[ind].ID}<br>{gff_obj[ind].locus_tag}<br>{gff_obj[ind].Name}<br>{gff_obj[ind].product}'
            ))
        
    fig = afns.update_fig_style(fig,
               title='Gene (selected from termination score graph)',
               ylabel = "Reads (per million)",
               xlabel = "Genomic position"
                       )
    return fig       


def graph_gene_region(df, feature_id, up_input, down_input,paths_to_valcnts_lst, 
                      path_annot_file,
                      annot_type='gbk'):
    """Plot the reads in region of the input feature with upstream and downstream genes.
    
    Arguments:
        df (pd.DataFrame): index should be the genomic position, each column 
        should be the value counts that we want mapped
        
    """
    #For generating annotation arrows etc
    start_offset, stop_offset, trimmed_df = afns.get_offsets(feature_id, up_input, 
                                        down_input, path_annot_file, annot_type=annot_type,
                                        contig_name=contig_name)
    l = []
    fig = go.Figure()
    #Actual read data
    for path in paths_to_valcnts_lst: 
        gp = gp.GeneProfile(path_annot_file, path)
        df, start, stop = gp.gene_hits(feature_id,
                             reads_per_mil=False,start_offset=start_offset,
                             stop_offset=stop_offset)
        fig.add_trace(go.Bar(
            x=df.index,
            y=df.loc[:,df.columns[0]],
            name = f'{contig_name}\n{feature_id}'
            ))
        l.append(df.loc[:,df.columns[0]].max())  
    max_val = max(l)
    #plot stop and start codon lines on graph
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
    #plot annoation arrows on graph using trimmed df
    for ind in trimmed_df.index:
        if trimmed_df.loc[ind,'strand'] == '-':
            stop = trimmed_df.loc[ind,'start']
            start = trimmed_df.loc[ind,'end']
        else:
            start = trimmed_df.loc[ind,'start']
            stop = trimmed_df.loc[ind,'end']   
        fig.add_trace(go.Scatter(
            x = pu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[0],
            y = pu.make_gene_arrow_coords(start, stop, max_val/10, -max_val/8)[1],
            fill="toself",
            mode = 'lines',
            line = {'width': 1.5,
                    'color':'rgb(200,200,200)',
                    },
            name=feature_id,
            showlegend=False,
            #text = f'{gff_obj[ind].ID}<br>{gff_obj[ind].locus_tag}<br>{gff_obj[ind].Name}<br>{gff_obj[ind].product}'
            ))
        
    fig = afns.update_fig_style(fig,
               title='Gene (selected from termination score graph)',
               ylabel = "Reads (per million)",
               xlabel = "Genomic position"
                       )
    return fig





