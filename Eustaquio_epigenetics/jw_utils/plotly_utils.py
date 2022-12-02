#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 28 09:50:33 2022

@author: jonwinkelman

"""
import numpy as np
from plotly import graph_objects as go
import plotly.offline as pyo


def plot_bar_with_outliers(series, name, end):
    """
    make a plotly 'histogram' with all values greater than 'end' binned into one bar
    
    parameters:
        series (pd.Series): input that will be plotted
        name (str): title of the graph
        end (int): all values greater than this int will be binned 
        
        return (go.Figure): 
    """
    start = int(series.min())
    binsize = 1
    # Making a histogram
    largest_value = series.max()
    if largest_value > end:
        hist = np.histogram(series, bins=list(range(start, end+binsize, binsize)) + [largest_value])
    else:
        hist = np.histogram(series, bins=list(range(start, end+binsize, binsize)) + [end+binsize])

    # Adding labels to the chart
    labels = []
    for i, j in zip(hist[1][0::1], hist[1][1::1]):
        if j <= end:
            labels.append('{}-{}'.format(i, j))
        else:
            labels.append('>{}'.format(i))

    # Plotting the graph
    trace = go.Bar(x=labels,
                   y=hist[0])

    return trace


def quick_scatter(x,y, mode = 'markers', plot=True):
    fig = go.Figure()
    fig.add_trace( go.Scatter(x = x,
                y = y,
                mode = mode
                         ))
    if plot:
        pyo.plot(fig)
    return fig

def quick_line(x,y, mode = 'lines', plot=True):
    fig = go.Figure()
    fig.add_trace( go.Scatter(x = x,
                y = y,
                mode = mode
                         ))
    if plot:
        pyo.plot(fig)
    return fig
    
    
def quick_bar(x,y, plot=True):
    fig = go.Figure()
    fig.add_trace( go.Bar(x = x,
                y = y,
                         ))
    if plot:
        pyo.plot(fig)
    return fig
   
def make_gene_arrow_coords(start_codon, stop_codon, height = 1,vert_shift = 0):
    """
    Return x and y arrays for arrow for gene
    Parameters:
    start_codon (int): genomic position of start codon
    stop_codon (int): genomic position of stop codon
    return (tuple of two lists)
    """
    if start_codon<stop_codon:
        strand = '+'
        br = stop_codon - (stop_codon - start_codon)/10
        print(br, start_codon, stop_codon)
        x = [start_codon,start_codon,br, br,stop_codon,br,br,start_codon,start_codon]
        y = [0, 0.5*height,0.5*height,height,0,-height,-0.5*height,-0.5*height,0]
    else: 
        strand = '-'
        br = stop_codon + (start_codon - stop_codon)/10 
        print(br, start_codon, stop_codon)
        x = [start_codon,start_codon,br, br,stop_codon,br,br,start_codon,start_codon]
        y = [0, 0.5*height,0.5*height,height,0,-height,-0.5*height,-0.5*height,0]
    if vert_shift:
        y = [ele+vert_shift for ele in y]
    return x,y
    
    