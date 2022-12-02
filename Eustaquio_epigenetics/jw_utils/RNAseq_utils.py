#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  1 14:03:09 2022

@author: jonwinkelman
"""
import pysam
import os
import pd
filepath = '/Users/jonwinkelman/Dropbox/Trestle_projects/Eustaquio_lab/Epigenetics_Project/RNAseq_bam_files/BAN-1.bam'
def get_bam_pairs(filename):
    samfile = pysam.AlignmentFile(filename, 'rb')  # BAM file reader.
    # Iterate through reads.
    read1 = None
    read2 = None
    
    pair_dict = {}
    for read in samfile:
        if read.flag==83:
            print(read.flag)
        if not read.is_paired or read.mate_is_unmapped or read.is_duplicate:
            continue
        if read.is_read2:
            read2 = read
        else:
            read1 = read
            continue
        if not read1 is None and not read2 is None and read1.query_name == read2.query_name:
             pair_dict[read1.query_name] = [read1, read2]
         #print(pair_dict[read1.query_name][0].flag ) 
    return pair_dict
             
min_length=0
max_length=1000
offset = 1 
pysam_arg = 'rb'
reverse_stranded = False       
fiveprime_end = True       
def basic_processing_fr_firststrand(min_length, max_length, filepath, offset = 16, 
                     pysam_arg = 'rb', fiveprime_end = True):
    pass
    """
    flags: 147 rna is on the '-' strand, is read2, fastq seq is rna seq
    Parameters:
    min_length (int): min length of the template that will be kept
    
    max_length (int): max length of the template that will be kept
        
    filepath (str): path to bam file aligned by bowtie2
    
    offset (int): number of bases back from template 3prime end, 1-based coordinates,
    
    
    Some notes: 
    r.seq is always + strand reference sequence, not necessarily the sequence of the read
    rna 3' refers to the template rna, not necessarily the read.
    flags 147 and 163 indicate that they are reads from the reverse sequencing reaction. IF , 
    thus the 3' of the template rna will be early in the read and of higher quality 
    
    *offset=15 selects the 15th residue back from the 3prime end of template rna, i.e. there
    are 14 residues that are ommited from the 3' end '
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
    new_path = f'./results/analysis/raw_filtered/{contig}_{new_filename}_{offset}_offset.txt'
    if os.path.isfile(new_path):
        raise Exception(f'the file {new_path} already exists')
    
    print(f'Creating file for {contig}')    
    samfile = samfilefull.fetch(contig=contig)  
    with open(new_path, 'w') as f:
        mapped = 0
        forward= 0
        proper_pair = 0
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
                        
                    # + strand genes      
                    elif r.flag == flags[1]:
                        f.write('\t'.join([ '+', r.seq, str(positions[0]), str(template_length) + '\n'  ]))
    source_name = filepath.split('/')[-1]                     
    #write log file           
    with open(f'./results/analysis/raw_filtered/{contig}_{new_filename}_{offset}.log', 'w') as f:
        f.write(f'an offset of {offset} was added into file\n')
        f.write(f'{source_name} contained {i} total reads\n')
        f.write(f'of {i} total reads, { (proper_pair/i)*100 }% were part of a proper pair\n')
        
    df = pd.DataFrame()
    data = [source_name, i, proper_pair, mapped]
    columns = ['source_name', 'offset', 'total_reads', 'proper_pair', 'mapped']
    for column, d in zip(columns, data):
        df[column] = [d]
df.to_csv(f'./results/analysis/raw_filtered/{new_filename}_{offset}.log.csv')
#return template_3            
