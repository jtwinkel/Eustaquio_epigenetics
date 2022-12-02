#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 15 13:12:48 2022

@author: jonwinkelman
"""

from jw_utils import file_utils as fu
from jw_utils import genome_utils as gu
from jw_utils import parse_gff as pgf
from jw_utils import ribo_profile3 as rp
from jw_utils import parse_gbk as pgb
import os
import pandas as pd
import numpy as np
import plotly.graph_objs as go



def check_for_proper_files(rp_object_dict, path_to_termination_dir, 
                           path_to_proteome, path_to_fasta_genome):
    #chceck for termination dir and files
    if not os.path.exists(path_to_termination_dir):
        os.makedirs(path_to_termination_dir)
    for exp in rp_object_dict.keys():
        filepath = os.path.join(path_to_termination_dir,f'{exp}_termscores.tsv')
        if not os.path.exists(filepath):
            df = rp_object_dict[exp].make_full_seq_df_term_scores(path_to_proteome,
                                path_to_fasta_genome, density_filt=0, 
                                density_display=True, reads_per_mil=True)
            with open(filepath, 'w') as f:
                cols = df.columns
                f.write(f'{df.index.name}\t{cols[0]}\t{cols[1]}\t{cols[2]}\t{cols[3]}\t{cols[4]}\n')
                for ind in df.index:
                    f.write(f'{ind}\t{df.loc[ind,cols[0]]}\t{df.loc[ind,cols[1]]}\t{df.loc[ind,cols[2]]}\t{df.loc[ind,cols[3]]}\t{df.loc[ind,cols[4]]}\n')
    #check for stop-start metagene dirs and files
    par_dirpath = './data/metagene_stops_starts'
    stop_dirpath = os.path.join(par_dirpath,'stops')
    start_dirpath = os.path.join(par_dirpath,'starts')
    if not os.path.exists(par_dirpath):
        os.makedirs(par_dirpath)
    if not os.path.exists(stop_dirpath):
        os.makedirs(stop_dirpath)
    if not os.path.exists(start_dirpath):
        os.makedirs(start_dirpath)
    
    for experiment in rp_object_dict.keys():
        stop_fpath_name = os.path.join(stop_dirpath, f'{experiment}_stop_norm.csv')
        start_fpath_name = os.path.join(start_dirpath, f'{experiment}_start_norm.csv')
        if not os.path.exists(stop_fpath_name):
            df = rp_object_dict[experiment].align_genes_stop(normalize=True)
            df = df.iloc[-100:,:]
            df.to_csv(stop_fpath_name)
        if not os.path.exists(start_fpath_name):    
            df = rp_object_dict[experiment].align_genes(normalize=True)
            df = df.iloc[:100,:]
            df.to_csv(start_fpath_name) 


def normalize_data(data):
    if (np.max(data) - np.min(data)) == 0:
        return np.full(len(data),0)
    else:
        return (data - np.min(data)) / (np.max(data) - np.min(data))
        

def linear_rescale(maximum, num_list):
    norm_data = normalize_data(num_list)
    norm_data = norm_data * maximum
    return norm_data


def make_term_score_fig(path_to_termination_dir, path_to_gff, 
                                    filt_value=0, clearance= -np.inf,
                                    return_type='fig'):   
    experimentalSetup_dict = simple_experimentalSetup_dict(path_to_termination_dir)
    cleared_genes = gu.get_genes_wo_neighbors(path_to_gff, clearance=clearance)
    fig = go.Figure()
    df_dict = {}
    for condition in experimentalSetup_dict.keys():
        df_lst = []
        i=0
        for replicate_path in experimentalSetup_dict[condition]:
            i+=1
            exp_name = expname_from_path(replicate_path)
            df = pd.read_csv(replicate_path, sep='\t', usecols=[0,1,2])
            df = df.set_index('protein ID')
            genes_to_keep = set(cleared_genes).intersection(set(df.index))
            df = df.loc[genes_to_keep,:]
            df.columns = [col+'_'+str(i) for col in df.columns]
            df_lst.append(df)
        df = pd.concat(df_lst,axis=1)
        dens_ave = (df.iloc[:,1] + df.iloc[:,3])/2
        df = df.loc[(dens_ave>=filt_value) , :]
        if return_type != 'fig':
            df_dict['-'.join(exp_name.split('-')[:2])] = df
        size = normalize_data(dens_ave)
        size=np.sqrt(size)
        size=np.sqrt(size)
        annot_df = pgf.make_simple_annot_df(path_to_gff).reset_index().set_index('protein_ID').loc[df.index,:]
        text = []
        for ind in annot_df.index:
            g = annot_df.loc[ind,'gene_ID']
            cn = annot_df.loc[ind,'common_name']
            p = annot_df.loc[ind,'product']
            text.append(f'{condition}<br>{ind}<br>{g}<br>{cn}<br>{p}')
        fig.add_trace(go.Scatter(
                        x = np.log2(df.iloc[:,0]),
                        y = np.log2(df.iloc[:,2]),
                        name = condition,
                        mode='markers',
                        marker = dict(
                                    size=size*40,
                                    #sizeref=2.*max(size)/(40.**2),
                                    sizemin=2
                                      ), 
                        text = text     
            ))
        fig.update_layout(
            title_text = 'Termination scores<br>(point size corresponds to read density)',
            paper_bgcolor='rgb(250,250,250)',   
            plot_bgcolor='rgb(253,253,253)',
            margin={'t':30, 'b':0, 'l':0, 'r':0}
            )
        fig.update_xaxes(
            title_text = "Replicate 1",
            title_font = {"size": 12},
            tickfont = {'size':8}
            )
        fig.update_yaxes(
            title_text = "Replicate 2",
            title_font = {"size": 12},
            tickfont = {'size':8}
            )
    if return_type == 'df':
        return df_dict
    elif return_type == 'all':
        return df_dict, fig 
    else: 
        return  fig   
            

def add_selections_to_feature_dropdown(path_to_gff):
    #get_all_feature names
    t = pgf.make_simple_annot_df(path_to_gff)
    
    return [{'label': f'{ind}, {t.loc[ind,"protein_ID"]}, name: {t.loc[ind,"common_name"]}, product: {t.loc[ind,"product"]}', 
             'value': t.loc[ind,"protein_ID"]} for ind in t.index]


def split_funct(path):
    'specific to path formatting, grabs the S1,S2 or S3'
    return '_'.join(path.split('/')[-1].split('-')[0:2])

def simple_experimentalSetup_dict(path_to_dir):
    """
    return a dictionary {condition SX: [paths_to_replicate1, 2, ...] }
    """
    val_counts_paths = fu.get_filepaths_in_dir(path_to_dir, end_to_exclude='.DS_Store')
    conditions = set([split_funct(path) for path in val_counts_paths])
    experimnent_setup_dict = {condition:[] for condition in conditions}
    for condition in conditions:
        for path in val_counts_paths:
            if condition == split_funct(path):
                experimnent_setup_dict[condition].append(path)
    return experimnent_setup_dict

def expname_from_path(path):
    'return experemental name part of filepath'
    return path.split('/')[-1].split('_')[0]




def filter_termScore_dfs(df_dict, high_cutoff=5, low_cutoff=0.2):
     'use output from make_termSeq_dfs_for_each_experiment()'
     high_df_dict = {condition:[] for condition in df_dict.keys()}
     low_df_dict = {condition:[] for condition in df_dict.keys()}
     for condition in df_dict.keys():
         for replicate_df in df_dict[condition]: # df_dict[condition] is a list of dfs
             high_filt = replicate_df['term scores'] > high_cutoff
             high_df_dict[condition].append(replicate_df[high_filt])
             low_filt = replicate_df['term scores'] < low_cutoff
             low_df_dict[condition].append(replicate_df[low_filt])
     return high_df_dict, low_df_dict


    

class rb_obj_store:
    'store the RiboProfile objects as attributes of this class that so '
    def __init__(self, list_of_valcount_paths, path_to_gff, feature_type='CDS'):
        
        self.path_to_gff = path_to_gff
        self.list_of_valcount_paths = list_of_valcount_paths
        self.feature_type = feature_type
        self.rp_object_dict = self.create_list_of_objects()
        self.keys = list(self.rp_object_dict.keys())
        
    
    def create_list_of_objects(self):
        rp_obj_dict = {}
        for path in self.list_of_valcount_paths:
            name = expname_from_path(path)
            rp_obj_dict[name] = rp.RiboProfile(self.path_to_gff, path, feature_type=self.feature_type)
        return rp_obj_dict


def get_trimmed_df(feature_ID, upstream, downstream, path_to_annot_file,
                     annot_type = 'gbk', contig_name=None):
    'get annotation df with the given # of up and downstream genes'
    if annot_type=='gbk':
        df = pgb.make_simple_annot_df(path_to_annot_file, contig_name, start_end=True).sort_values('start')
    else:
        df = pgf.make_simple_annot_df(path_to_annot_file, start_end=True).sort_values('start')
    
    #df = df.set_index('protein_ID')
    index = df.index.get_loc(feature_ID)
    if index-upstream<=0:
        upstream = index
    max_index = df.index.get_loc(df.index.max())
    if index+downstream>=max_index:
        downstream = max_index-index
    trimmed_df = df.iloc[index-upstream:(index+downstream)+1,:]
    return trimmed_df


def get_offsets(feature_id, up_input, down_input, path_to_annot_file,
                 annot_type = 'gbk', contig_name=None):
    trimmed_df = get_trimmed_df(feature_id, up_input, down_input,
                         path_to_annot_file, annot_type=annot_type, 
                         contig_name=contig_name)
    """
    Return the nt positions of start and end of region, as well as the trimmed df.
    
    Arguments:
    feature_id (str): ID of feature, needs to match the index of the trimmed_df. If using 
                      genbank annotation files, I think should use locus tags for feature IDs. 
    up_input (int): number of genes to calculate offsets for  upstream of the feature
    down_input (int): number of genes to calculate offsets for downstream of the feature
    path_to_annot_file (str): 
    annot_type (str): default='gbk', can also be gff, not tested here
    contig_name (str): default=None, name of contig, chromosome, plasmid...
    
    """

    s = trimmed_df.loc[trimmed_df.index.min(),'start'] # start of first upstream gene
    e = trimmed_df.loc[trimmed_df.index.max(),'end'] # end of last downstream gene
    feature_start = trimmed_df.loc[feature_id,'start']
    feature_end = trimmed_df.loc[feature_id,'end']
    feature_phase = trimmed_df.loc[feature_id,'strand']
    if feature_phase == '+':
        start_offset = feature_start-s
        stop_offset = e-feature_end
    if feature_phase == '-':
        stop_offset = feature_start-s
        start_offset = e-feature_end
    if start_offset==0:start_offset=50
    if stop_offset==0:stop_offset=50
    return start_offset, stop_offset, trimmed_df
    
    
    
    
def density_filter(filt_level):
    'get read_density for each gene and return set of genes that have higher density than input'
    exp_dict = simple_experimentalSetup_dict('./data/termination_files')
    filt_genes_dict = {}
    for condition in exp_dict.keys():  
        for replicate in exp_dict[condition]:
            name = expname_from_path(replicate)
            df = pd.read_csv(replicate, sep='\t', usecols=[0,1,2,3])
            filt = df.loc[:,df.columns[2]] >= filt_level
            filt_genes_dict[name] = list(df.loc[filt,df.columns[0]])  
    return filt_genes_dict
    
    
def stop_start_traces(sums_ser, line_dict, name):
    trace = go.Scatter(
                x = list(sums_ser.keys()),
                y=sums_ser.values,
                mode = 'lines',
                line = line_dict,
                name = name,
                text = name
                )
    return trace
    
    
def plot_stops_starts_metagene(path_to_gff, stops_or_start='starts', density_filt =False, 
                               filt_level=0, neighbor_clearance=0,
                               return_type='fig'):
    exp_dict = simple_experimentalSetup_dict(f'./data/metagene_stops_starts/{stops_or_start}')
    
    neigh_filt_features = gu.get_genes_wo_neighbors(path_to_gff, clearance=neighbor_clearance)
    fig = go.Figure()
    df_dict= {}
    if density_filt:
        filt_genes_dict = density_filter(filt_level)
    for condition in exp_dict.keys():  
        for replicate in exp_dict[condition]:
            df = pd.read_csv(replicate)
            df = df.set_index(df.columns[0])
            print(f'genes before neigh_filt: {df.shape[1]}')
            df = df.loc[:,neigh_filt_features]
            print(f'genes after neigh_filt: {df.shape[1]}')
            name = expname_from_path(replicate)
            if density_filt:
                filt_genes = filt_genes_dict[name]
                filt_genes = set(filt_genes).intersection(set(neigh_filt_features))
                df = df.loc[:,filt_genes]
                print(f'genes after density_filt: {df.shape[1]}')
            print(f'final # genes: {df.shape[1]}')
            if return_type !='fig':
                df_dict[name] = df
            else:
                df_means = df.mean(axis=1)
                line_dict = {'width':3}
                trace = stop_start_traces(df_means, line_dict, name=name)
                fig.add_trace(trace)
    if return_type == 'fig': return fig      
    else: return df_dict   


def update_fig_style(fig, title=None, ylabel=None,xlabel=None):
    fig.update_layout(
        title_text = title,
        paper_bgcolor='rgb(250,250,250)',   
        plot_bgcolor='rgb(253,253,253)',
        margin={'t':30, 'b':0, 'l':0, 'r':0}
        )
    fig.update_xaxes(
        title_text = xlabel,
        title_font = {"size": 12},
        tickfont = {'size':8}
        )
    fig.update_yaxes(
        title_text = ylabel,
        title_font = {"size": 12},
        tickfont = {'size':8}
        )
    return fig