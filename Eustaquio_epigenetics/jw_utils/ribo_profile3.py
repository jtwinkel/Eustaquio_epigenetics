
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 17 14:18:21 2022

@author: jonwinkelman
"""
import pandas as pd
import bisect
from plotly import offline as pyo
from plotly import graph_objects as go
import numpy as np
import os
from jw_utils import parse_gff as pgf
from jw_utils import parse_fasta as pf
from jw_utils import genome_utils as gu

class RiboProfile:
    
    def __init__(self, path_to_gff, path_to_value_counts, feature_type = 'CDS'):
        self.path_to_gff = path_to_gff
        self.feature_type = feature_type
        self.genome_size = self.get_genome_descriptors()['genome_size']
        self.gene_objects = pgf.make_seq_object_dict(path_to_gff, feature_type) 
        self.genes = list(self.gene_objects.keys())
        self.path_to_value_counts = path_to_value_counts
        self.value_counts = self.read_value_counts()

        self.plus_strand_genes = [gene for gene in self.genes if self.gene_objects[gene].strand == '+']
        self.minus_strand_genes = [gene for gene in self.genes if self.gene_objects[gene].strand == '-']
        self.total_mapped_reads = self.total_mapped_reads()
        self.gff_geneObject_dict = pgf.make_seq_object_dict(
                                            self.path_to_gff, 
                                            feature_type=feature_type
                                            )
                                            
    def hit_density(self, feature_id, start_offset=0, stop_offset=0, reads_per_mil=False):
        'return the (#total reads in gene region )/ (length of the gene + offsets)'
        g_obj = self.gff_geneObject_dict[feature_id]
        length = (g_obj.end - g_obj.start) + (start_offset + stop_offset)
        hits = self.gene_hits(feature_id, start_offset=start_offset, stop_offset=stop_offset, 
                        reads_per_mil=reads_per_mil)[0]
        return hits.sum()/length
    
    
    def get_reference_sequence(self, path_to_fasta_reference):
        return pf.get_seq_dict(path_to_fasta_reference)
        
        
    def get_reference_sequence(self, path_to_fasta_reference):
        return pf.get_seq_dict(path_to_fasta_reference)
        
        
        
    def termination_scores(self, gene_list = None, add = .01, reads_per_mil = False,
                           body_length_normalization=True, density_filt=0,df_return=False):
        '''
        Parameters:
            density_filt (int): number of hits/nt (if reads_per_mil = False)
                if reads_per_mil = True, then this is 
                (total reads in gene in reads_per_mil)/genelength
        Body = Body of ORF minus 30 nucleotides at 5’ and 3’ end
        Tail = Last 9 nucleotides of ORF, corresponding to last two sense codons + STOP
        term_score =tail_hits/body_hits
        Returns: 
            term_scores: (dict)
            dictionary of {feature ID : termination score}
            

        '''

        term_scores = {}
        flags = {}
        if not gene_list:
            gene_list = self.gff_geneObject_dict.keys()
        for protein in gene_list:
            hits_df, start_codon, stop_codon = self.gene_hits_start_zero(
                                                        protein, 
                                                        start_offset=0, 
                                                        stop_offset=0, 
                                                        reads_per_mil=reads_per_mil
                                                        )
            gene_length = stop_codon - start_codon
            density = hits_df.loc[:,hits_df.columns[0]].sum()/hits_df.shape[0]
            if density >= density_filt:
                if gene_length > 79:
                    body_hits_df = hits_df.loc[start_codon+30:stop_codon-30,hits_df.columns[0]]
                    tail_hits_df = hits_df.loc[stop_codon-8:stop_codon,hits_df.columns[0]]    
                    body_hits = body_hits_df.sum()
                    tail_hits = tail_hits_df.sum()
                    if body_hits == 0:
                        body_hits = add
                    if tail_hits == 0:
                        tail_hits = add
                    if body_length_normalization:
                        body_hits = body_hits/body_hits_df.shape[0]
                        tail_hits = tail_hits/tail_hits_df.shape[0]
                        term_score = tail_hits/body_hits
                        b = '{:.2f}'.format(body_hits)
                        t = '{:.2f}'.format(tail_hits)
                        ts ='{:.2f}'.format(term_score)
                        flags[protein] =f'<br>body rds/mil/nt: {b}<br>tail rds/mil/nt: {t}<br>term_score: {ts}<br>'
                    else:
                        term_score = tail_hits/body_hits
                        b = '{:.2f}'.format(body_hits)
                        t = '{:.2f}'.format(tail_hits)
                        ts ='{:.2f}'.format(term_score)
                        flags[protein] =f'<br>body rds/mil: {b}<br>tail rds/mil: {t}<br>term_score: {ts}<br>'
                    if df_return == True:
                        term_scores[protein] = [tail_hits/body_hits, density]
                    else:
                        term_scores[protein] = tail_hits/body_hits
                        
        if df_return == True:
            df = pd.DataFrame(term_scores).transpose()
            if reads_per_mil:
                columns = ['term scores','read density, r/mil']
            else:
                columns = ['term scores','read density']
            df.columns = columns
            return df
        else:
            return term_scores
                    
                    
        


    
    def get_genome_descriptors(self):
        desc_dict = {}
        with open(self.path_to_gff, 'r') as f:
            for i,line in enumerate(f):
                if line.startswith('#!genome-build-accession'):
                    desc_dict['assembly_accession'] = line.split(':')[-1].strip()
                if line.startswith('#!annotation-date'):
                    desc_dict['annotation_date'] = ' '.join(line.split(' ')[-2:]).strip()
                
                if line.startswith('##sequence-region'):
                    desc_dict['sequence_region'] = ' '.join(line.split(' ')[-3:]).strip()
                if line[0] != '#' and line[0] != ' ':
                    anot_lst = line.split('\t')
                    if anot_lst[2] == 'region':
                        desc_dict['genome_size'] = int(anot_lst[4])
                        break
        return desc_dict
  
        
    def read_value_counts(self):
        'return counts of 3 prime ends at each genomoic position'
        df = pd.read_csv(self.path_to_value_counts)
        df.columns = ['genomic_position', "3'_(+)", "3''_(-)"]
        df = df.fillna(0)
        df = df.set_index('genomic_position')
        return df#df.reindex(list(range(1,self.genome_size+1)),fill_value=0)

    def total_mapped_reads(self):
        'return total # of hits on the reference genome'
        plus_counts = self.value_counts.iloc[:,0].sum()
        minus_counts = self.value_counts.iloc[:,1].sum()
        return plus_counts + minus_counts
    
    def NormalizeData(self, data):
        if (np.max(data) - np.min(data)) == 0:
            return np.full(len(data),0)
        else:
            return (data - np.min(data)) / (np.max(data) - np.min(data))
            
    
    def gene_hits(self, gene_name, start_offset=50, stop_offset=50, normalize = False, reads_per_mil = False):
        df = self.value_counts
        gene_obj = self.gene_objects[gene_name]
        if gene_obj.strand == '+' :
            start_codon = gene_obj.start
            stop_codon = gene_obj.end
            start_index = bisect.bisect_left(df.index,(start_codon - start_offset))
            end_index = bisect.bisect_right(df.index, (stop_codon + stop_offset)) 
            df2 = df.iloc[start_index:end_index,0]

        elif gene_obj.strand == '-':
            start_codon = gene_obj.end
            stop_codon = gene_obj.start
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
    
    
    def gene_hits_with_nulls(self, gene_name, df=None, start_offset=50, stop_offset=50, normalize = False, reads_per_mil = False):
        gene_obj = self.gene_objects[gene_name]
        
        if type(df) != pd.DataFrame:
            df=self.value_counts
            df = df.reindex(list(range(1,self.genome_size+1)),fill_value=0)
        
        if gene_obj.strand == '+' :
            start_codon = gene_obj.start
            stop_codon = gene_obj.end
            beg =gene_obj.start - start_offset
            end = gene_obj.end + stop_offset
            df2 = df.loc[beg:end,df.columns[0]]

        elif gene_obj.strand == '-':
            start_codon = gene_obj.end
            stop_codon = gene_obj.start
            beg = gene_obj.start - stop_offset
            end =gene_obj.end + start_offset
            df2 = df.loc[beg:end,df.columns[1]]
    
        
        if reads_per_mil:
            mil_reads = self.total_mapped_reads/1000000
            df2 = df2.div(mil_reads)
        if normalize:
            if (np.max(df2) - np.min(df2)) != 0:
                df2 = self.NormalizeData(df2)

        return pd.DataFrame(df2), start_codon, stop_codon
        
        
        
    
    def gene_hits_start_zero(self, gene, df=None, start_offset=50, 
                            stop_offset=50, normalize=False, reads_per_mil=False):
        if type(df) != pd.DataFrame:
            df=self.value_counts
            df = df.reindex(list(range(1,self.genome_size+1)),fill_value=0)
        gene_obj = self.gene_objects[gene]
        df2, start_codon, stop_codon = self.gene_hits_with_nulls(gene, 
                                                 df = df,
                                                 start_offset=start_offset, 
                                                 stop_offset=stop_offset, normalize=normalize,
                                                 reads_per_mil=reads_per_mil)
        
        gene_length = gene_obj.end - gene_obj.start
        if gene_obj.strand == '+':
            df2.index = df2.index- gene_obj.start
            start_codon = start_codon - gene_obj.start
            stop_codon = stop_codon - gene_obj.start
            
        if gene_obj.strand == '-':
            df2.index = df2.index - gene_obj.end #
            start_codon = start_codon - gene_obj.end
            stop_codon = start_codon + gene_length
            df2.index = df2.index*-1
            df2 = df2.sort_index()
        return df2, start_codon, stop_codon
    
    def genes_hit_sums(self, start_offset=0, stop_offset=0, reads_per_mil=False, density=False):
        'return sum of hits within the gene'
        gene_sums = {}
        for gene in self.genes:
            df = self.gene_hits(gene, start_offset=start_offset, stop_offset=stop_offset, 
                                reads_per_mil=reads_per_mil)[0]
            gene_sums[gene] = df.loc[:,df.columns[0]].sum(axis=0)
            
        return pd.DataFrame.from_dict(gene_sums,orient ='index', columns = ['sum of hits'])
    
    
    def gene_hits_stop_zero(self, gene, df=None, start_offset=50, 
                            stop_offset=50, normalize=False, reads_per_mil=False):
        if type(df) != pd.DataFrame:
            df=self.value_counts
            df = df.reindex(list(range(1,self.genome_size+1)),fill_value=0)
        gene_obj = self.gene_objects[gene]
        df2, start_codon, stop_codon = self.gene_hits_with_nulls(gene, 
                                                 df = df,
                                                 start_offset=start_offset, 
                                                 stop_offset=stop_offset, normalize=normalize,
                                                 reads_per_mil=reads_per_mil)
        
        if gene_obj.strand == '+':
            df2.index = df2.index- gene_obj.end
            start_codon = start_codon - gene_obj.end
            stop_codon = stop_codon - gene_obj.end
        if gene_obj.strand == '-':
            df2.index = df2.index- gene_obj.start
            start_codon = (start_codon - gene_obj.start)*-1
            stop_codon = stop_codon -  gene_obj.start
            df2.index = df2.index*-1
            df2 = df2.sort_index()
        return df2, start_codon, stop_codon
    
    
        
    def plot_gene(self, gene, start_aligned=False, stop_aligned=False, 
                  normalize = False, plot = True, reads_per_mil=False):
        gene_obj = self.gene_objects[gene]
        if start_aligned:
            df2, start_codon, stop_codon = self.gene_hits_start_zero(gene, normalize = normalize, reads_per_mil=reads_per_mil)
        elif stop_aligned:
            df2, start_codon, stop_codon = self.gene_hits_stop_zero(gene, normalize = normalize, reads_per_mil=reads_per_mil)
        else:
            df2, start_codon, stop_codon = self.gene_hits(gene, normalize = normalize, reads_per_mil=reads_per_mil)
        if plot:
            y = [0, 0.25,0.25,0.5,0,-0.5,-0.25,-0.25,0]
            max_val = df2.loc[:,df2.columns[0]].max()
# =============================================================================
#             if gene_obj.strand == '-' and not stop_align_flipped and not start_align_flipped:  
#                 br = stop_codon + abs(stop_codon - start_codon)/10
#             else:
# =============================================================================
            br = stop_codon - abs(stop_codon - start_codon)/10
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
            fig.add_trace(go.Scatter(
                x = [start_codon, start_codon, br, br, stop_codon,br,br,start_codon,start_codon], 
                y = [ele*(max_val/10)-max_val/15 for ele in y],
                line=dict(color="rgb(100,100,100)"),
                name = 'strand: '+gene_obj.strand))
            pyo.plot(fig)
        return df2, start_codon, stop_codon



    def align_genes(self,genelist=None, normalize=False, reads_per_mil=False,
                    start_offset=50,stop_offset=50):
        'return a df with all genes aligned at start (default) or stop codons'
        valcount_df = self.value_counts
        valcount_df = valcount_df.reindex(list(range(1,self.genome_size+1)),fill_value=0)
        gene_plots = []
        if not genelist:
            genelist = self.genes
        for gene in genelist:
            #df, start, stop = self.plot_gene(gene, start_aligned=True, plot=False, reads_per_mil=reads_per_mil)
            df, start, stop = self.gene_hits_start_zero(gene,
                                        df = valcount_df,
                                        start_offset=start_offset,
                                        stop_offset=stop_offset, 
                                        normalize=normalize,
                                        reads_per_mil=reads_per_mil)
            df.columns = [gene]   
            gene_plots.append(df)
        return pd.concat(gene_plots, axis=1, join = 'outer').fillna(0)
    
    
    
    def align_genes_stop(self,genelist=None, normalize=False, reads_per_mil=False,
                         start_offset=50,stop_offset=50):
        'return a df with all genes aligned at start (default) or stop codons'
        valcount_df = self.value_counts
        valcount_df = valcount_df.reindex(list(range(1,self.genome_size+1)),fill_value=0)
        gene_plots = []
        if not genelist:
            genelist = self.genes
        for gene in genelist:
            df, start, stop = self.gene_hits_stop_zero(gene, 
                                                   df =valcount_df,
                                                   start_offset=start_offset,
                                                   stop_offset=stop_offset, 
                                                   normalize=normalize,
                                                   reads_per_mil=reads_per_mil)
            df.columns = [gene]   
            gene_plots.append(df)
        return pd.concat(gene_plots, axis=1, join = 'outer').fillna(0)
            
         
        
    
    
    def plot_start_aligned_genes(self, genelist=None, normalize=False, min_reads=0,
                                 return_type='trace_only', reads_per_mil=False, 
                                 return_full_df=False, start_offset=50,
                                 stop_offset=50):
        'plot sum etc. of all genes aligned at start'
        if not genelist:
            genelist=self.genes
        
        df = self.align_genes(genelist=genelist, 
                           start_offset=start_offset, stop_offset=stop_offset, 
                           normalize=normalize, reads_per_mil=reads_per_mil)
        lst_of_normcols =[]
        if min_reads > 0:                  
            for col in df.columns:
                if df.loc[:,col].sum()>min_reads:
                    lst_of_normcols.append(col)
            df=df.loc[:,lst_of_normcols]
                    
        df = df.mean(axis=1)

        trace = go.Scatter(
            x=df.index[:100],  
            y=df[:100],
            mode = 'lines',
            name = self.path_to_value_counts.split('/')[-1].split('_')[0],
            #line=dict(color="rgb(0,255,0)")
            )
        if return_full_df:
            return_type = ''
            return full_df
        if return_type == 'trace_only':
            return trace
        if return_type  =='df':
            fig = go.Figure()
            fig.add_trace(trace)
            pyo.plot(fig)
            return df
    
    
    def plot_stop_aligned_genes(self, genelist=None, normalize = False, min_reads = 0, 
    return_type='trace_only', reads_per_mil=False, return_full_df=False):
        'plot sum etc. of all genes aligned at stop'
        
        if normalize:
            df = self.align_genes_stop(genelist = genelist, reads_per_mil=reads_per_mil)
            lst_of_normcols =[]
            for col in df.columns:
                total_reads = df.loc[:,col].sum()
                if total_reads>min_reads:
                    lst_of_normcols.append(self.NormalizeData(df.loc[:,col]))
            full_df = pd.concat(lst_of_normcols, axis=1)
            df = full_df.mean(axis=1)
        else:
            full_df = self.align_genes_stop(genelist = genelist, reads_per_mil=reads_per_mil)
            df = full_df.mean(axis=1)
            
        trace = go.Scatter(
            x =df.index[-100:], 
            y=df[-100:],
            mode = 'lines',
            name = self.path_to_value_counts.split('/')[-1].split('_')[0],
            #line=dict(color="rgb(0,255,0)")
            )
        if return_full_df:
            return_type = ''
            return full_df
        if return_type == 'trace_only':
            return trace
        if return_type  =='df':
            fig = go.Figure()
            fig.add_trace(trace)
            pyo.plot(fig)
            return df
        
        
        
        
    def make_full_seq_df_term_scores(self, path_to_proteome,
                                    path_fasta_genome, density_filt=0.25, 
                                    density_display=False, reads_per_mil=False):
        
        '''
        return termination scores and the associated nt and AA seqs and ids
        parameters:
            density_filter (int): do not include any genes with fewer than 
            this average number of reads per nt (total reads/gene length)
        
        '''
        
        if self.feature_type !='CDS':
            raise Exception(f'feature type for ribo_profile2 object is {self.feature_type}\n'
                            'it must be set to "CDS"')
        term_scores = self.termination_scores(density_filt=density_filt)
        if density_display:
            term_scores_df = self.termination_scores(density_filt=density_filt,
                                                        df_return=True, reads_per_mil=reads_per_mil)
            term_scores_df.index.rename('protein ID', inplace=True)
            
            
        else:
            term_scores = self.termination_scores(density_filt=density_filt)
            term_scores_df = pd.DataFrame(
                    {'protein ID':term_scores.keys(),'term scores':term_scores.values()}
                    ).set_index('protein ID')
        gu_object = gu.ProteomeUtils(path_to_proteome, self.path_to_gff)
        df_DNA_seqs = gu_object.get_DNA_seq_for_proteins(path_fasta_genome,term_scores.keys())
        df_AA_seqs = gu_object.get_get_AA_seq_for_proteins(term_scores.keys())
        
        df_t = term_scores_df.merge(df_DNA_seqs, on = 'protein ID')
        return df_t.merge(df_AA_seqs, on = 'protein ID')

        
    









    
  