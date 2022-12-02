#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 27 14:39:32 2022

@author: jonwinkelman
"""

from jw_utils import parse_gff as pgf
import pandas as pd
from jw_utils import parse_fasta as pf
import numpy as np


def get_genes_wo_neighbors(path_to_gff, clearance=-5000, return_df=False):
    ''''
    get genes that have > clearance separation from neighbors on same strand
    
    parameters:
        path_to_gff (str): path to the gff3 genome annotation file
        clearance (int): number of nt that must exist beteween gene and 
            neighbor before start and after stop codon (on the same dna strand)
        return_df (bool): if False, returns list of genes that meet clearance
            criteria. if True, data-rich df is returned
        
        returns (list or df): see return_df parameter   
    '''
    
    gff_dict = pgf.make_seq_object_dict(path_to_gff,feature_type='CDS')
    strands = ['+', '-']
    df_lst = []
    for strand in strands:
        starts = []
        ends = []
        genes = []
        phase = []
        gene_df = pd.DataFrame()
        for gene in gff_dict.keys():
            gene = gff_dict[gene]
            if gene.strand == strand:
                starts.append(gene.start)
                ends.append(gene.end)
                genes.append(gene.ID)
                phase.append(gene.strand)
        gene_df['feature_ID'] = genes    
        gene_df['start'] = starts 
        gene_df['end'] = ends 
        gene_df['phase'] = phase 
        gene_df = gene_df.sort_values('start', ascending=True).reset_index(drop=True)
        gene_df['start_clearance'] = gene_df['start']         - gene_df['end'].shift()
        gene_df['end_clearance'] = gene_df['start'].shift(-1) - gene_df['end']
        gene_df = gene_df.fillna(np.inf)
        filt = (gene_df.loc[:,'start_clearance']>clearance) & (gene_df.loc[:,'end_clearance'] >clearance)
        gene_df = gene_df.loc[filt,:]
        df_lst.append(gene_df)
    if return_df:
        return  pd.concat(df_lst).set_index('feature_ID')
    else:
        return list(df_lst[0]['feature_ID']) + list(df_lst[1]['feature_ID'])



def get_genes_wo_neighbors_df(path_to_gff, clearance=-5000):
    ''''
    get genes that have > clearance separation from neighbors on same strand
    
    parameters:
        path_to_gff (str): path to the gff3 genome annotation file
        clearance (int): number of nt that must exist beteween gene and 
            neighbor before start and after stop codon (on the same dna strand)
            
        returns (df): Data-rich df is returned
    '''
    return get_genes_wo_neighbors(path_to_gff, clearance=clearance, return_df=True)
    


def rev_comp(seq):
    ''' 
    return the reverse complement of a dna sequences

    This is not optimized for large sequences, but code is simple
    
    Parameters:
    
    seq (str): can be a string of upper or lower case DNA bases: (a,t,c,g)
    
    
    '''
    base_pair_dict = {
        'A':'T',
        'T':'A',
        'C':'G',
        'G':'C',
        'a':'t',
        't':'a',
        'c':'g',
        'g':'c',
        }
    rev_comp_seq = ''
    for base in seq:
        rev_comp_seq = base_pair_dict[base] + rev_comp_seq
    return rev_comp_seq 


def get_gene_seq(path_to_fasta_genome):
    fasta_dict = pf.get_seq_dict(path_to_fasta_genome) 
    
    
class GenomeUtils:
    def __init__(self, path_to_fasta_genome, path_to_gff):
         
        self.path_to_fasta_genome = path_to_fasta_genome
        self.genome_seq_dict = pf.get_seq_dict(self.path_to_fasta_genome) 
        self.path_to_gff = path_to_gff
        self.fasta_dict = pf.get_seq_dict(self.path_to_fasta_genome)
        self.gff_cdsObject_dict = pgf.make_seq_object_dict(
                                            self.path_to_gff, 
                                            feature_type='CDS'
                                            ) 
        self.gff_geneObject_dict = pgf.make_seq_object_dict(
                                            self.path_to_gff, 
                                            feature_type='gene'
                                            ) 
        self.chromosome_names = list(self.genome_seq_dict.keys())

        
        
    def get_gene_dna_seq(self, gene):
        polarity = self.gff_geneObject_dict[gene].strand
        gene_start = self.gff_geneObject_dict[gene].start
        gene_end = self.gff_geneObject_dict[gene].end
        gene_seq = self.genome_seq_dict[self.chromosome_names[0]][(gene_start-1):gene_end]
        if polarity == '-':
            start_codon = gene_end
            stop_codon = gene_start
            gene_seq = rev_comp(gene_seq)
        return gene_seq
        
        

        
class ProteomeUtils:
    def __init__(self, path_to_proteome, path_to_gff):  
        self.path_to_proteome = path_to_proteome
        self.path_to_gff = path_to_gff
        self.protein_seq_dict =  pf.get_seq_dict(self.path_to_proteome)
     
    
    def get_AA_seq_for_protein(self, protein):
        return self.protein_seq_dict[protein]
        
    def get_get_AA_seq_for_proteins(self, protein_list):
        AA_seqs = [self.protein_seq_dict.get(protein.replace('cds-','')) for protein in protein_list]
        return pd.DataFrame({'protein ID':protein_list, 'AA_seq':AA_seqs}).set_index('protein ID')
            
            
    def get_DNA_seq_for_proteins(self, path_fasta_genome, protein_list):
        '''
        returns the DNA sequences for the input proteins
        parameters]
    
        '''
        genome_obj = GenomeUtils(path_fasta_genome, self.path_to_gff)
        prot2gene_dict = pgf.make_prot2gene_dict(self.path_to_gff)
        seqs = []
        genes = []
        proteins = []
        for protein in protein_list:
            gene = prot2gene_dict.get(protein)
            if genome_obj.gff_geneObject_dict.get(gene):
                seq = genome_obj.get_gene_dna_seq(gene)
            else:
                seq = ''
            seqs.append(seq)
            genes.append(gene)
            proteins.append(protein)
        df = pd.DataFrame({'protein ID':proteins,'gene ID':genes, 'sequence':seqs})
        return df.set_index('protein ID')

        
    
            
            
            
    
