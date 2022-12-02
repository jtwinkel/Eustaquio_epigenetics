#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 08:44:54 2022

@author: jonwinkelman
"""
import pysam
import pandas as pd

def pileup_on_each_contig(bam_filepath, write_file=False):
    """
    return a df with number of reads at each position in each chromosome
    parameters:
        bam_filepath (str): should point to a sorted? bam file
        
    return (DataFrame):
    df structure: 
        index: genomic postion (one-based numbering)
        columns: each column represents one chromosome, plasmid, etc
        
        index starts at 0, the first nucleotide, and ends at the last nucleotide
        of the longest chromosome, plasmid, etc
        
        if there are no reads a genomic position within the length of the contig,
        then the reads will be equal to 0. If the genomic postion is out of the
        range of contig, then read value at that position will be equal to np.nan

    """
    bam_file = pysam.AlignmentFile(bam_filepath, 'rb')
    lst_of_dfs = []
    if len(bam_file.references) >1:
        for contig in bam_file.references:
            pileup_genome = {}
            for col in bam_file.pileup(contig = contig):
                if pileup_genome.get(contig):
                    pileup_genome[contig].update({col.reference_pos:col.nsegments})
                else:
                    pileup_genome[contig] = {col.reference_pos:col.nsegments}
            chromosome_length = bam_file.get_reference_length(contig)
            df = pd.DataFrame(pileup_genome)
            df = df.reindex(index=range(chromosome_length), fill_value=0)
            lst_of_dfs.append(df)
    else:
        contig = bam_file.references[0]
        pileup_genome = {}
        for col in bam_file.pileup():
            if pileup_genome.get(contig):
                pileup_genome[contig].update({col.reference_pos:col.nsegments})
            else:
                pileup_genome[contig] = {col.reference_pos:col.nsegments}
        chromosome_length = bam_file.get_reference_length(contig)
        df = pd.DataFrame(pileup_genome)
        df = df.reindex(index=range(chromosome_length), fill_value=0)
        lst_of_dfs.append(df)
        
    df = pd.concat(lst_of_dfs, axis=1, join='outer')
    df.index = df.index +1 #change to 1-based numbering
    df.index.names = ['position']
    if write_file:
        name = bam_filepath.split("/")[-1]
        df.to_csv(f'./{name}_pilup.csv')
    return df