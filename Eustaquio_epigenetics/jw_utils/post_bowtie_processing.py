#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 11 12:17:02 2022

@author: jonwinkelman
"""
import pysam
import pandas as pd
import numpy as np
import os
from jw_utils import file_utils


def basic_processing(min_length, max_length, filepath, offset = 16, path=None):
    """
    
    Parameters:
    min_length (int): min length of the template that will be kept
    
    max_length (int): max length of the template that will be kept
        
    filepath (str): path to bam file aligned by bowtie2
    
    offset (int): number of bases back from template 3prime end, 1-based coordinates,
    
    
    Some notes: 
    r.seq is always + strand reference sequence, not necessarily the sequence of the read
    rna 3' refers to the template rna, not necessarily the read.
    flags 147 and 163 indicate that they are reads from the reverse sequencing reaction, 
    thus the 3' of the template rna will be early in the read and of higher quality 
    
    *offset=15 selects the 15th residue back from the 3prime end of template rna, i.e. there
    are 14 residues that are ommited from the 3' end '
    """
    
    if offset == 0:
        raise Exception('there is no zeroth residue to select, this is 1-based numbering')
    samfile = pysam.AlignmentFile(filepath, 'rb')
    new_filename = filepath.split('/')[-1].split('.')[0] +  '_filtered'
    if not os.path.isdir('./results/analysis'):
        os.makedirs('./results/analysis')
    if not os.path.isdir('./results/analysis/raw_filtered'):
        os.makedirs('./results/analysis/raw_filtered')
    if not path:
        new_path = f'./results/analysis/raw_filtered/{new_filename}_{offset}_offset.txt'
    else:
        new_path = path
    if os.path.isfile(new_path):
        raise Exception(f'the file {new_path} already exists')
        
    with open(new_path, 'w') as f:
        filtered = 0
        revcomps =0
        mapped = 0
        forward= 0
        proper_pair = 0
        raw_lengths = []
        f.write('polarity' + '\t' + 'sequence' + '\t' + 'rna_3prime'+ '\t' + 'temp_len' + '\n')
        for i, r in enumerate(samfile):
            if r.is_proper_pair:
                proper_pair +=1
                if r.is_mapped:
                    mapped +=1
                    template_length = abs(r.tlen)
                    raw_lengths.append(template_length)
                    if template_length >= min_length & template_length <= max_length:
                        filtered +=1
                        positions = r.get_reference_positions()
                        positions = [p+1 for p in positions] # convert to 1-based numbering
                        # + strand genes
                        if r.flag == 147:  #r.seq is reverse complement of actual read   
                            f.write('\t'.join([ '+', r.seq, str(positions[-(offset)]), str(template_length) + '\n'  ]))
                            revcomps+=1
                        # - strand genes      
                        elif r.flag == 163:
                            f.write('\t'.join([ '-', r.seq, str(positions[offset-1]), str(template_length) + '\n'  ]))
                            revcomps+=1 
                        elif r.flag == 99:
                            forward+=1
                        elif r.flag == 83:
                            forward+=1
    source_name = filepath.split('/')[-1]                     
    #write lengths_to_file
    with open(f'./results/analysis/raw_filtered/{new_filename}_{offset}.lengths', 'w') as f:
        for length in raw_lengths:
            f.write(str(length) + '\n')
    
    #write log file           
    with open(f'./results/analysis/raw_filtered/{new_filename}_{offset}.log', 'w') as f:
        f.write(f'an offset of {offset} was added into file\n')
        f.write(f'{source_name} contained {i} total reads\n')
        f.write(f'of {i} total reads, {filtered} ({filtered/i}) were part of aligned pairs longer than {min_length} and shorter than {max_length}\n')
        f.write(f'of {filtered} filtered reads, {mapped} ({mapped/filtered}) were mapped to the genome\n')
        f.write(f'of {filtered} filtered & mapped reads, {revcomps} ({revcomps/filtered}) were from the reverse squencing reaction\n')
        f.write(f'of {filtered} filtered & mapped reads, {forward} ({forward/filtered}) were from the forward squencing reaction\n')
        #f.write(f'of {mapped} filtered & mapped reads, {proper_pair} ({proper_pair/mapped}%) were part of a proper pair.  ')
    return i, filtered, mapped, revcomps            


def basic_processing_chromosomes(read_obj_list, chr_name, min_length, max_length, filepath, offset = 16, path=None):
    """
    
    Parameters:
    read_obj_list (list): list of samfile read objects from a single chromosome
    chr_name (str): name of the chromosome, plasmid etc...
    min_length (int): min length of the template that will be kept
    max_length (int): max length of the template that will be kept   
    filepath (str): path to bam file aligned by bowtie2
    offset (int): number of bases back from template 3prime end, 1-based coordinates,
    
    
    Some notes: 
    r.seq is always + strand reference sequence, not necessarily the sequence of the read
    rna 3' refers to the template rna, not necessarily the read.
    flags 147 and 163 indicate that they are reads from the reverse sequencing reaction, 
    thus the 3' of the template rna will be early in the read and of higher quality 
    
    *offset=15 selects the 15th residue back from the 3prime end of template rna, i.e. there
    are 14 residues that are ommited from the 3' end '. If one is selected, it is simply the 
    last nt on the 3' end of the RNA, the first nt in read of the reverse reaction.
    """
    
    if offset == 0:
        raise Exception('there is no zeroth residue to select, this is 1-based numbering')
    
    new_filename = filepath.split('/')[-1].split('.')[0] +  '_filtered_' + chr_name
    if not os.path.isdir('./results/analysis'):
        os.makedirs('./results/analysis')
    if not os.path.isdir('./results/analysis/raw_filtered'):
        os.makedirs('./results/analysis/raw_filtered')
    if not path:
        new_path = f'./results/analysis/raw_filtered/{new_filename}_{offset}_offset.txt'
    else:
        new_path = path
    if os.path.isfile(new_path):
        raise Exception(f'the file {new_path} already exists')
        
    with open(new_path, 'w') as f:
        filtered = 0
        revcomps =0
        mapped = 0
        forward= 0
        proper_pair = 0
        raw_lengths = []
        f.write('polarity' + '\t' + 'sequence' + '\t' + 'rna_3prime'+ '\t' + 'temp_len' + '\n')
        for i, r in enumerate(read_obj_list):
            if r.is_proper_pair:
                proper_pair +=1
                if r.is_mapped:
                    mapped +=1
                    template_length = abs(r.tlen)
                    raw_lengths.append(template_length)
                    if template_length >= min_length & template_length <= max_length:
                        filtered +=1
                        positions = r.get_reference_positions()
                        positions = [p+1 for p in positions] # convert to 1-based numbering
                        # + strand genes
                        if r.flag == 147:  #r.seq is reverse complement of actual read   
                            f.write('\t'.join([ '+', r.seq, str(positions[-(offset)]), str(template_length) + '\n'  ]))
                            revcomps+=1
                        # - strand genes      
                        elif r.flag == 163:
                            f.write('\t'.join([ '-', r.seq, str(positions[offset-1]), str(template_length) + '\n'  ]))
                            revcomps+=1 
                        elif r.flag == 99:
                            forward+=1
                        elif r.flag == 83:
                            forward+=1
    source_name = filepath.split('/')[-1]                     
    #write lengths_to_file
    with open(f'./results/analysis/raw_filtered/{new_filename}_{offset}.lengths', 'w') as f:
        for length in raw_lengths:
            f.write(str(length) + '\n')
    
    #write log file           
    with open(f'./results/analysis/raw_filtered/{new_filename}_{offset}.log', 'w') as f:
        f.write(f'an offset of {offset} was added into file\n')
        f.write(f'{source_name} contained {i} total reads\n')
        f.write(f'of {i} total reads, {filtered} ({filtered/i}) were part of aligned pairs longer than {min_length} and shorter than {max_length}\n')
        f.write(f'of {filtered} filtered reads, {mapped} ({mapped/filtered}) were mapped to the genome\n')
        f.write(f'of {filtered} filtered & mapped reads, {revcomps} ({revcomps/filtered}) were from the reverse squencing reaction\n')
        f.write(f'of {filtered} filtered & mapped reads, {forward} ({forward/filtered}) were from the forward squencing reaction\n')
        #f.write(f'of {mapped} filtered & mapped reads, {proper_pair} ({proper_pair/mapped}%) were part of a proper pair.  ')
    return i, filtered, mapped, revcomps 



def sorted_counts_df(filepath, offset = 16, write_to_file = True):
    '''
    return sorted df with number of rna 3 prime (maybe offset) ends at each genomic position
    
    Parameters:
    
    filepath (str): path to filtered offset file output by basic_processing function
    
    returns:
        
    dataframe 
    also writes dataframe to csv file
    ''' 
    df = pd.read_csv(filepath, sep = '\t', usecols = [0,2])
    print(f'DF has been read. shape is {df.shape}')
    df_r = df[(df.iloc[:,0] == '-')]
    df_f = df[(df.iloc[:,0] == '+')]
    print(f'DF has been split. shapes are: df_f= {df_f.shape},  df_r = {df_r.shape}')
    df_f = df_f.iloc[:,1].value_counts(sort=False).sort_index()
    df_r = df_r.iloc[:,1].value_counts(sort=False).sort_index()
    print(f'DF has been value counted. shapes are: df_f= {df_f.shape},  df_r = {df_r.shape}')
    df = pd.merge(df_f, df_r, right_index=True, left_index=True,how='outer')
    print(f'DF has been value merged. shapes is: {df.shape}')
    col1 = f"3'_{offset}_offset_(+)"
    col2 = f"3'_{offset}_offset_(-)"
    df.columns = [col1, col2]
    new_filename = filepath.split('/')[-1]
    parent_dir = './results/analysis/value_counts'
    if not os.path.isdir(parent_dir):
        os.makedirs(parent_dir)
    new_filepath = f'./results/analysis/value_counts/{new_filename}_val_cnts.csv'
    if os.path.isfile(new_filepath):
        raise Exception(f'the file {new_filepath} already exists')
    if write_to_file:
        df.to_csv(new_filepath)
    return df








