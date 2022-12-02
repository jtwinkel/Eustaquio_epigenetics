#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  8 20:34:50 2022

@author: jonwinkelman
"""
from Bio import SeqIO

out_path = './reference.gff'

def gb_to_gff(path_to_genbank, out_path):
    """Convert genbank annotation into gff3 annotation.
    
    A work in progress - many fields still need to be added
    """
    gb_obj = SeqIO.parse(path_to_genbank, "genbank")
    with open(out_path, 'w') as fp:
        fp.write('##gff-version 3\n')
        for contig in gb_obj:  
            fp.write(f'{contig.id}\tUnknown\tregion\t{str(1)}\t{str(len(contig.seq))}\t.\t.\t.\tID={contig.id};') #not sure what to put in the 'strand' category
            fp.write(f'taxon={contig.annotations["taxonomy"]};')
            for f in contig.features:
                l = f.location
                if l.strand == -1: strand='-'
                else: strand ='+'
                fq = f.qualifiers
                if not f.type=='source':
                    fp.write(f'{contig.id}\tUnknown\t{f.type}\t{int(l.start)}\t{int(l.end)}\t.\t{strand}\t.\t')
                if f.type == 'CDS':
                    pid = fq.get("protein_id")
                    if pid: pid = pid[0]
                    fp.write(f'ID=cds-{pid};Parent=gene-{fq["locus_tag"][0]};product={fq["product"][0]}\n')
                if f.type == 'gene':
                    gid = fq.get("locus_tag")
                    if gid: gid = gid[0]
                    fp.write(f'ID=cds-{gid};locus_tag={gid}\n')
                if f.type == 'rRNA':
                    gid = fq.get("locus_tag")
                    if gid: gid = gid[0]
                    fp.write(f'ID=cds-{gid};Parent=gene-{fq["locus_tag"][0]};product={fq["product"][0]}\n')
                if f.type == 'ncRNA':
                    gid = fq.get("locus_tag")
                    if gid: gid = gid[0]
                    fp.write(f'ID=cds-{gid};Parent=gene-{fq["locus_tag"][0]};product={fq["product"][0]}\n')
                if f.type == 'tRNA':
                    gid = fq.get("locus_tag")
                    fp.write(f'ID=cds-{gid};Parent=gene-{fq["locus_tag"][0]};product={fq["product"][0]}\n')

