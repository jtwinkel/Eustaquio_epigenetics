#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:03:09 2022

@author: jonwinkelman
"""
import pysam
import os
import pandas as pd
from jw_utils import file_utils as fu


path_to_dir = '/Users/jonwinkelman/Dropbox/Trestle_projects/Eustaquio_lab/new_alignments/Bowtie_aligned'
paths = fu.get_filepaths_in_dir(path_to_dir, end_to_exclude='.DS_Store')




def get_bam_pairs(fpath_to_bam, contig=None):
    """Return two dicts where."""
    samfile = pysam.AlignmentFile(fpath_to_bam, 'rb')  # BAM file reader.
    if contig:
        samfile = samfile.fetch(contig=contig)
    d_read1 = {}
    d_read2 = {}
    for read in samfile:
        if read.is_read1:
            d_read1[read.query_name] = read
        elif read.is_read2:
            d_read2[read.query_name] = read
    df1 = pd.DataFrame.from_dict(d_read1, orient='index')
    df2 = pd.DataFrame.from_dict(d_read2, orient='index')
    df1.columns = ['r1_objects']
    df2.columns = ['r2_objects']
    df1 = pd.merge(df1, df2, right_index=True, left_index=True)
    del df2
    return df1


def get_flag_dist(path_to_bam):
    flag_vals = range(0,1000)
    flag_d = {flag:0 for flag in flag_vals}
    samfile = pysam.AlignmentFile(path_to_bam, 'rb')  # BAM file reader.
    for r in samfile:
        flag_d[r.flag] +=1
    df = pd.DataFrame(flag_d.keys(), flag_d.values())
    df.columns = ['counts']
    df = df.loc[df['counts'] >0,:]
    return df
            
            

             
min_length=0
max_length=1000
offset = 1 
pysam_arg = 'rb'
reverse_stranded = False       
fiveprime_end = True       

def basic_processing_fr_firststrand(filepath, 
                     pysam_arg = 'rb', fiveprime_end = True):
    pass
    """
    flags: 147 rna is on the '-' strand, is read2, fastq seq is rna seq
    Parameters:
        
    filepath (str): path to bam file aligned by bowtie2
    
    
    
    Some notes: 
    r.seq is always + strand reference sequence, not necessarily the sequence of the read
    rna 3' refers to the template rna, not necessarily the read.
    flags 147 and 163 indicate that they are reads from the reverse sequencing reaction. IF , 
    thus the 3' of the template rna will be early in the read and of higher quality 
    
    """
    if fiveprime_end:
        flags = [147,163]
    else:
        flags = [99,83]
    if offset == 0:
        raise Exception('there is no zeroth residue to select, this is 1-based numbering')
    
    samfilefull = pysam.AlignmentFile(filepath, pysam_arg)
    new_filename = filepath.split('/')[-1].split('.')[0] +  '_filtered'
    if not os.path.isdir('./results/analysis'):
        os.makedirs('./results/analysis')
    if not os.path.isdir('./results/analysis/raw_filtered'):
        os.makedirs('./results/analysis/raw_filtered')
        
    for contig in samfilefull.references:
        new_path = f'./results/analysis/raw_filtered/{contig}_{new_filename}.txt'
        if os.path.isfile(new_path):
            raise Exception(f'the file {new_path} already exists')  
        print(f'Creating file for {contig}')    
        samfile = samfilefull.fetch(contig=contig)  
        with open(new_path, 'w') as f:
            mapped = 0
            proper_pair = 0
            total_hits = 0
            f.write('polarity' + '\t' + 'sequence' + '\t' + 'rna_5prime'+ '\t' + 'temp_len' + '\n')
            for i, r in enumerate(samfile):
                template_length = abs(r.tlen)
                if r.is_proper_pair:
                    proper_pair +=1               
                    if r.is_mapped:
                        mapped +=1
                        positions = r.get_reference_positions()
                        positions = [p+1 for p in positions] # convert to 1-based numbering
                        # - strand genes
                        if r.flag == flags[0]:  #r.seq is reverse complement of actual read   
                            f.write('\t'.join([ '-', r.seq, str(positions[-1]), str(template_length) + '\n'  ]))
                            total_hits +=1
                            
                        # + strand genes      
                        elif r.flag == flags[1]:
                            f.write('\t'.join([ '+', r.seq, str(positions[0]), str(template_length) + '\n'  ]))
                            total_hits +=1
        source_name = filepath.split('/')[-1]                     
        #write log file           
        with open(f'./results/analysis/raw_filtered/{contig}_{new_filename}.log', 'w') as f:
            f.write(f'{source_name} contained {i} total reads\n')
            f.write(f'of {i} total reads, { (proper_pair/i)*100 }% were part of a proper pair\n')
            f.write(f'there were { (proper_pair/2)} complete pair sets\n')
            f.write(f'of { (proper_pair/2)} pair sets, {total_hits} 5prime ends were mapped\n')
            
        df = pd.DataFrame()
        data = [source_name, i, proper_pair/2, total_hits]
        columns = ['source_name','total_reads', 'proper_pairs', 'total_hits']
        df_dict = {}
        for column, d in zip(columns, data):
            df[column] = [d]
            df.to_csv(f'./results/analysis/raw_filtered/{new_filename}.log.csv')
            df_dict[contig] = df
    return df  



def sorted_counts_df(filepath, write_to_file = True):
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
    col1 = "5'_(+)"
    col2 = "5'_(-)"
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
