#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 09:29:33 2022

@author: jonwinkelman
"""
from operator import length_hint
import pandas as pd
import bisect
from plotly import offline as pyo
from plotly import graph_objects as go
import numpy as np
import os
from scipy import ndimage as ndi
from jw_utils import app_functions as afns
from jw_utils import parse_gff as pgf
from jw_utils import parse_fasta as pf
from jw_utils import genome_utils as gu
from jw_utils import parse_gbk as pgb
from jw_utils import plotly_utils as pu
from jw_utils import file_utils as fu
from Bio import SeqIO



class GeneProfile:
    """Class for mapping reads to features in a contig using a genbank file."""
    
    def __init__(self, path_to_annot_file, path_to_value_counts=None,
                 val_counts_df = None, contig_name=None):
        """Contig object for mapping reads to features with a genbank file.
        
        Arguments:
            path_to_annot_file (str): Path the the genbank annotaion file
            
            path_to_value_counts (str): path to csv contianing genomic position 
            in the index and the counts of 5' ends of reads at each genomic 
            position on the top strand: "5'_(+)_sum" and bottom strand: "5'_(-)_sum"
            
            val_counts_df (DataFrame): df contianing genomic position 
            in the index and the counts of 5' ends of reads at each genomic 
            position on the top strand: "5'_(+)_sum" and bottom strand: "5'_(-)_sum"
        
        """


        if path_to_value_counts:
            self.value_counts = self.read_value_counts()
            self.contig_name = path_to_value_counts.split('/')[-1].split('_')[0]
        else:
            self.value_counts = val_counts_df
        if type(val_counts_df) == pd.core.frame.DataFrame:
            if not contig_name:
                raise Exception('Need value for the argument "contig_name"')
            else:
                self.contig_name = contig_name
        self.path_to_value_counts = path_to_value_counts
        self.path_to_annot_file = path_to_annot_file
        self.geneObject_dict = self.build_genbank_dict()[self.contig_name]
        self.genes = list(self.geneObject_dict.keys())
        self.plus_strand_genes = [gene for gene in self.genes if self.geneObject_dict[gene].strand == 1]
        self.minus_strand_genes = [gene for gene in self.genes if self.geneObject_dict[gene].strand == -1]
        self.total_mapped_reads = self.total_mapped_reads()
        
        
        self.genome_size = self.get_genome_descriptors()['genome_size']
    
    def get_genome_descriptors(self):
        """Return a dict with diff genome descriptors."""
        description_dict = {}
        gb_obj = SeqIO.parse(self.path_to_annot_file, "genbank")
        for chromosome in gb_obj:
            if chromosome.name == self.contig_name:
                description_dict['genome_size'] = len(chromosome.seq)
        return description_dict
    
    
    
    def build_genbank_dict(self):
        """Build dict {chromosomeID:{featureID:gbk_feature_object}} from a genbank annotation file."""
        chromosome_dict = {}
        gb_obj = SeqIO.parse(self.path_to_annot_file, "genbank")
        for chromosome in gb_obj:  
            for feature in chromosome.features:
                if feature.type == 'CDS':
                    if chromosome_dict.get(chromosome.name): 
                        chromosome_dict[chromosome.name].update({feature.qualifiers['locus_tag'][0]:feature})
                    else:
                        chromosome_dict[chromosome.name] = {feature.qualifiers['locus_tag'][0]:feature}
        return chromosome_dict



    def read_value_counts(self):
        """Return counts of 3 prime ends at each genomoic position."""
        df = pd.read_csv(self.path_to_value_counts)
        df.columns = ['genomic_position', "5'_(+)", "5''_(-)"]
        df = df.fillna(0)
        df = df.set_index('genomic_position')
        return df#df.reindex(list(range(1,self.genome_size+1)),fill_value=0)



    def total_mapped_reads(self, *args, **kwargs):
        """Return total number of reads present in value counts df."""
        df = self.value_counts
        return df[df.columns[0]].sum() + df[df.columns[1]].sum()


    def NormalizeData(self, data):
        """Scale data between 0 and 1, unless all values are 0, then return all zeros."""
        if (np.max(data) - np.min(data)) == 0:
            return np.full(len(data),0)
        else:
            return (data - np.min(data)) / (np.max(data) - np.min(data))


    def gene_hits(self, gene_name, start_offset=100, stop_offset=100, reads_per_mil=False, normalize=False):
        """Map reads from value counts file to a gene region."""
        df = self.value_counts
        gene_obj = self.geneObject_dict[gene_name]
        if gene_obj.strand == 1 :
            start_codon = gene_obj.location.start
            stop_codon = gene_obj.location.end
            start_index = bisect.bisect_left(df.index,(start_codon - start_offset))
            end_index = bisect.bisect_right(df.index, (stop_codon + stop_offset)) 
            df2 = df.iloc[start_index:end_index,0]

        elif gene_obj.strand == -1:
            start_codon = gene_obj.location.end
            stop_codon = gene_obj.location.start
            start_ind = bisect.bisect_left(df.index,(stop_codon - stop_offset))
            end_ind = bisect.bisect_right(df.index, (start_codon + start_offset)) 
            df2 = df.iloc[start_ind:end_ind,1]
    
            
        if reads_per_mil:
            mil_reads = self.total_mapped_reads/1000000
            df2 = df2.div(mil_reads)
        if normalize:
            if (np.max(df2) - np.min(df2)) != 0:
                df2 = self.NormalizeData(df2)

        return pd.DataFrame(df2), start_codon, stop_codon
    
    
    def plot_gene(self, gene,  start_offset=100, stop_offset=100, 
                  normalize = False, plot = True, reads_per_mil=False,
                  vert_shift=-1):
        """
        Return df, start and stop codon pos, plot gene hits with the gene arrow.
        
        Arguments:
        gene (str)
        normalize (bool)
        plot (bool)
        reads_per_mil (bool)
        """
        gene_obj = self.geneObject_dict[gene]
        if gene_obj.strand == -1: strand = '-'
        else: strand = '+'
        df2, start_codon, stop_codon = self.gene_hits(gene, start_offset=start_offset, 
                                            stop_offset=stop_offset,normalize = normalize, 
                                            reads_per_mil=reads_per_mil)
        if plot:
            max_val = df2.loc[:,df2.columns[0]].max()
            fig = go.Figure()
            fig.add_trace(go.Scatter(
                x=df2.index,
                y=df2.loc[:,df2.columns[0]],
                name = gene))
            fig.add_trace(go.Scatter(
                x=[start_codon,start_codon], 
                y=[0,max_val],
                line=dict(color="rgb(0,255,0)"),
                name = 'Start'))
            fig.add_trace(go.Scatter(
                x=[stop_codon,stop_codon], 
                y=[0,max_val],
                line=dict(color="rgb(255,0,0)"),
                name = 'Stop'))
            x,y = pu.make_gene_arrow_coords(start_codon, stop_codon, height = 1,vert_shift = vert_shift)
            fig.add_trace(go.Scatter(
                x = x, 
                y = y,
                line=dict(color="rgb(100,100,100)"),
                name = 'strand: '+strand))
            pyo.plot(fig)
            return fig
        return df2, start_codon, stop_codon
    
    
    
    def plot_gene_region(self, feature_id, up_input, down_input,
                          annot_type='gbk', convolve=False, kernal=None):
        """Plot the reads in region of the input feature with upstream and downstream genes.
        
        Arguments
        - feature_id (str): gene name as formatted in annotation file
        - up_input (int) number of upstream genes to plot
        -down_input (int) number of downstream genes to plot
        
        Return
        - plotly Figure object
        """
        fig = go.Figure()
        #if up_input or down_input:     
        #gff_obj=pgf.make_seq_object_dict(path_to_gff, feature_type='CDS')
        
        start_offset, stop_offset, trimmed_df = afns.get_offsets(feature_id, up_input, 
                                            down_input, self.path_to_annot_file,
                                            annot_type=annot_type,
                                            contig_name=self.contig_name)
        l = []
            
        df, start, stop = self.gene_hits(feature_id,
                                 reads_per_mil=True,start_offset=start_offset,
                                 stop_offset=stop_offset)
        fig.add_trace(go.Bar(
            x=df.index,
            y=df.loc[:,df.columns[0]],
            name = f'{self.contig_name}\n{feature_id}'
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







class TssFinder(GeneProfile):
    """Class to find transcription start sites."""
    
    def get_tss_windows(self, window_size=100):
        """Return dictionary with genome coordinates for windows likely to house the tss.
        Arguments:
            window_size (int): number of nucleotides before start codon to consider for tss
        Return: dict structure: {feature_id:beginning,end} beginning and end is always relative
                        to the top strand of the genome, e.g the beginning still is a 
                        smaller number than end even on '-' strand.
        """
        windows = {}
        #get tss window for genes on the plus strand
        for i, gene in enumerate(self.plus_strand_genes):
            start_codon = int(self.geneObject_dict[gene].location.start)
            window_start = start_codon - window_size
            windows[gene] = window_start, start_codon
            #determine if upstream gene overlaps into tss window, if so then set window coordinates to None
            if i>0:
                upstream_gene = self.plus_strand_genes[i-1]
                if self.geneObject_dict[upstream_gene].location.end > window_start:
                    windows[gene] = None
        #get tss window for genes on the minus strand
        for i, gene in enumerate(self.minus_strand_genes):
            start_codon = int(self.geneObject_dict[gene].location.end)+1
            window_end = start_codon + window_size
            windows[gene] =  start_codon, window_end
            if i < len(self.minus_strand_genes)-1:
                downstream_gene = self.minus_strand_genes[i+1]
                if self.geneObject_dict[downstream_gene].location.start < window_end:
                    windows[gene] = None
        return windows
    
    
    def get_window_read_density(self, window_size=100):
        """Return the number of reads/per/nt in the tss window.
        
        Arguments
        - window_size (int): length of the window to look for tss
            upstream of the start codon
            
        Return
        - dictionary with each gene as a key and reads/nt as the value
        """
        windows = self.get_tss_windows(window_size=window_size)
        for gene in windows:
            if windows[gene]:
                li = bisect.bisect_left(self.value_counts.index, windows[gene][0])
                ri = bisect.bisect_right(self.value_counts.index, windows[gene][1])
                if self.geneObject_dict[gene].strand==1:
                    df = self.value_counts.iloc[li:ri,0]
                else:
                    df = self.value_counts.iloc[li:ri,1]
                windows[gene] = df.sum()/window_size
        return windows
    
    
    
    def plot_window_region(self,feature_id, up_input, down_input, window_size=100):
        windows = self.get_tss_windows(window_size=window_size)
        feat_window_coords = windows.get(feature_id)
        fig = self.plot_gene_region(feature_id, up_input, down_input)
        fig.add_vrect(x0 =feat_window_coords[0], x1 = feat_window_coords[1],
                      line_width=0, fillcolor="red", opacity=0.2)
        pyo.plot(fig)
        return fig
 
    
 
    def fill_df(self, max_val=None):
        df=self.value_counts
        if not max_val:
            max_val = self.genome_size
        return df.reindex(range(1,max_val),fill_value=0)
 
 
    def get_region_hits(self, window_size=100):
        """Return a dict of window dfs for gene 5' tss window with the reads at each position"""
        windows = self.get_tss_windows(window_size=window_size)
        df = self.fill_df()
        gene_windows = {}
        for feature in windows.keys():
            feat_window_coords = windows[feature]
            if feat_window_coords:
                #li = bisect.bisect_left(self.value_counts.index, feat_window_coords[0])
                #ri = bisect.bisect_right(self.value_counts.index, feat_window_coords[1])
                if self.geneObject_dict[feature].strand==1:
                    gene_windows[feature] = df.loc[feat_window_coords[0]:feat_window_coords[1],df.columns[0]]
                    #gene_windows[feature] = self.value_counts.iloc[li:ri,0]
                else:
                    gene_windows[feature] = df.loc[feat_window_coords[0]:feat_window_coords[1],df.columns[1]]
                    #gene_windows[feature] = self.value_counts.iloc[li:ri,1]
        return gene_windows
    
    def smooth_df(self, kernal_type, kernal_length):
        """Return a df smoothed with a chosen kernal
        
        Arguments:
            kernal_type (str): gauss or mean
            length (int): lenght of the kernal
        """
        for col in self.value_counts.columns:
            #kernal
            pass
            #s = ndi.correlate(df[col], )
    

def get_zero_centered_range(length):
    """Return np array of input length that is centered around 0."""
    xi =np.arange(length)
    x0 = length//2
    return xi-x0    

    
def gauss_1d_kernal(length, sigma=1, plot=False):
    """Return a 1 dimensional gaussian kernal."""
    if length%2==0:
        raise Exception('The kernal length must be an odd number')
    x = get_zero_centered_range(length)
    sigma = sigma
    gaussian_kernal = (1/(np.sqrt(2*np.pi)*sigma) * 
                      np.exp(-(x**2) / 2*sigma**2) )
    if plot:
        pu.quick_scatter(x, gaussian_kernal)
    return gaussian_kernal


def mean_diff_kernal(length, plot=False):
    """Return a 1 dimensional kernal that mean smooths and finds edge."""
    if length%2==0:
        raise Exception('The kernal length must be an odd number')
    
    mean_kernal = np.full(length,  1/length)
    diff_kernal = [-1,0,1]
    mean_diff_kernal = np.convolve(mean_kernal, diff_kernal, mode = 'full')
    x = get_zero_centered_range(len(mean_diff_kernal))
    if plot:
        pu.quick_scatter(x, mean_diff_kernal)
    return mean_diff_kernal


def gausian_diff_kernal(length, sigma, plot=False):
    """Return a 1 dimensional kernal that gaussian smooths and finds edge."""
    if length%2==0:
        raise Exception('The kernal length must be an odd number')
    x = get_zero_centered_range(length)  
    gaussian_kernal = gauss_1d_kernal(length, sigma=sigma, plot=False)
    diff_kernal = [-1,0,1]
    gausian_diff_kernal = np.convolve(gaussian_kernal, diff_kernal, mode = 'full')
    x = get_zero_centered_range(len(gausian_diff_kernal))
    if plot:
        pu.quick_scatter(x, gausian_diff_kernal)
    return gausian_diff_kernal



def make_fig_dict(kernal_type='gauss_diff', length=7, sigma=1)
    df_dict = {}
    fig_dict = {}
    kernal = [1,0,-1]
    if kernal_type = 'gauss_diff'
        kernal = gpro.gausian_diff_kernal(length=length, sigma=sigma, plot=True)
    for gene in tss_win_dict.keys():
        df_dict[gene] = ndi.correlate(tss_win_dict[gene],kernal, mode='reflect')
        fig = pu.quick_line(x = tss_win_dict[gene].index, y = df_dict[gene], plot=False)
        af = pu.quick_bar(x= tss_win_dict[gene].index, y =tss_win_dict[gene], plot=False)['data'][0]
        fig.add_trace(af)
        fig_dict[gene] = fig
    return fig_dict

