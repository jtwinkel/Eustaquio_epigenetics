#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 21 20:38:09 2022

@author: jonwinkelman
"""
import pysam
import post_bowtie_processing as pbp
filepath = '/Users/jonwinkelman/Dropbox/Trestle_projects/Eustaquio_lab/Epigenetics Project/RNAseq_bam_files/BAN-1.bam'
# =============================================================================
# samfile = pysam.AlignmentFile(filepath, 'rb')
# chroms = ['BF000000.1','BF000000.2','BF000000.3','BF000000.4','BF000000.5']
# chrom_dict = {chrom:[] for chrom in chroms}
# for r in samfile:
#     if r.reference_name:
#         chrom_dict[r.reference_name].append(r)
# =============================================================================



pbp.basic_processing_chromosomes(chrom_dict['BF000000.5'], 'BF000000.5', 15,500,filepath, offset=1)