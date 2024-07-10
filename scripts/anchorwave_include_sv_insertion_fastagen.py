#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:03:42 2022

@author: nathan
"""

import sys
import pysam
import re
from collections import defaultdict


"""
Usage:
    
    python anchorwaave_include_sv_insertion_fastagen.py \
        anchorwave_tsv_file.tsv 
        ref_genotype \
        query_genotype \
        ref_genome.fasta \
        query_genome.fasta


"""

def extract_sv_plus_adjacent_seqs(ins_geno, non_ins_geno, 
                                  insertion_fasta_file, absent_insertion_fasta_file, 
                                  ins_chromosome_num, rowID, ins_start, ins_end,
                                  non_ins_coord, non_ins_chromosome, left_aln_region, right_aln_region):
    """
    Extracts SV + 300bp flanking alignable regions in insertion genotype and 300 bp flanking alignable region in opposite genotype where insertion is absent
    """
    # left_fasta = list()
    # right_fasta = list()
    """First, extract SV+300bpflanking alignable regions in insertion genotype
    """
    if left_aln_region != 'NA' and int(float(left_aln_region)) >= 300 and right_aln_region != 'NA' and int(float(right_aln_region)) >= 300:

        """
        Extract all seq (left 300bp + SV + right 300bp)
        """
        start = int(float(ins_start) - 300)
        end = int(float(ins_end) + 300)
        left_sv_right = pysam.faidx(
            insertion_fasta_file, ins_chromosome_num+':'+str(start)+'-'+str(end))
        ins_fasta_re_change_fasta_id = re.sub(
            r'\>\S+', '>'+ins_geno+'_ins_chr'+ins_chromosome_num+':'+str(start)+'-'+str(end)+'_'+non_ins_genotype+':'+non_ins_coord+'_'+rowID, left_sv_right)
        # return left_sv_right_rowID_fasta

        # left_seq_chunk_end = int(float(ref_start) - 1)

        # left_chunk_fasta = pysam.faidx(insertion_fasta_file, 'chr'+chromosome_num+':'+str(start)+'-'+str(end))
        # left_rowID_fasta = re.sub(r'\>\S+', '>B73_chr'+chromosome_num+'_ins_'+rowID, left_chunk_fasta)

        """
        Extact 300 bp flanking alignable region at insertion location in genome WITHOUT insertion
        """

        non_ins_start = int(float(non_ins_coord) - 300)
        non_ins_end = int(float(non_ins_coord) + 300)
        non_ins_fasta = pysam.faidx(absent_insertion_fasta_file, non_ins_chromosome+':'+str(non_ins_start)+'-'+str(non_ins_end))
        non_ins_fasta_re_change_fasta_id = re.sub(r'\>\S+', '>'+non_ins_geno+'_chr'+non_ins_chromosome+':'+str(
            non_ins_start)+'-'+str(non_ins_end)+'_'+rowID, non_ins_fasta)
        return ins_fasta_re_change_fasta_id, non_ins_fasta_re_change_fasta_id
    else:
        ins_fasta_re_change_fasta_id = 'NA'
        non_ins_fasta_re_change_fasta_id = 'NA'
        return ins_fasta_re_change_fasta_id, non_ins_fasta_re_change_fasta_id

# def write_chunk_fasta(genotype, chromosome_num, left_fasta, right_fasta):
#     with open(genotype+'_insertion_chr'+chromosome_num+'_chunks.fasta', 'w') as of:
#         for i in left_fasta:
#             if i != 'NA':
#                 of.write(i)
#         for j in right_fasta:
#             if j != 'NA':
#                 of.write(j)


#################################################################################################################################################################################################################
#################################################################################################################################################################################################################
#################################################################################################################################################################################################################

# "/home/nathan/Downloads/chr1_B73_Oh7B_StructuralVariant_Summary.tsv"
chr_tsv_file = sys.argv[1]
flanking_sv_dict = defaultdict(list)
counter = 0
output_fasta = list()
non_ins = list()
left_fasta_query = list()
right_fasta_query = list()
chromosome_num_ref = ''
chromosome_num_query = ''
genotype_ref = ''
genotype_query = ''
with open(chr_tsv_file, 'r') as fh:
    next(fh)
    for line in fh:
        # counter += 1
        line2 = line.strip()
        line3 = line2.split()
        # rowID = line3[12]
        Block_Type = line3[6]
        chromosome = line3[0]
        if Block_Type.split('_in')[-1] == sys.argv[2]:

            # chromosome_num_ref = line3[0]
            ins_genotype = sys.argv[2]
            non_ins_genotype = sys.argv[3]
            ins_fasta = sys.argv[4]
            non_ins_fasta = sys.argv[5]

            
            ins_chromosome_num = line3[0]
            rowID = line3[12]
            # Block_Type = current_line_list[6]
            ins_start = line3[1]
            ins_end = line3[2]
            left_aln_region = line3[14]
            right_aln_region = line3[17]
            non_ins_coord = line3[4]
            non_ins_chromosome = line3[3]

            
            ins, non_ins = extract_sv_plus_adjacent_seqs(ins_genotype, non_ins_genotype, ins_fasta, non_ins_fasta, ins_chromosome_num, rowID, ins_start, ins_end, non_ins_coord, non_ins_chromosome, left_aln_region, right_aln_region)
            output_fasta.extend([ins, non_ins])

        else:

            ins_genotype = sys.argv[3]
            non_ins_genotype = sys.argv[2]
            ins_fasta = sys.argv[5]
            non_ins_fasta = sys.argv[4]            
            ins_chromosome_num = line3[3]
            rowID = line3[12]
            # Block_Type = current_line_list[6]
            ins_start = line3[4]
            ins_end = line3[5]
            left_aln_region = line3[14]
            right_aln_region = line3[17]
            non_ins_coord = line3[1]
            non_ins_chromosome = line3[0]            
            ins, non_ins = extract_sv_plus_adjacent_seqs(ins_genotype, non_ins_genotype, ins_fasta, non_ins_fasta, ins_chromosome_num, rowID, ins_start, ins_end, non_ins_coord, non_ins_chromosome, left_aln_region, right_aln_region)

            output_fasta.extend([ins, non_ins])

with open(sys.argv[2] + '_' + sys.argv[3] + '_chr'+chromosome+'_pseudoref.fasta', 'w') as of:
    for i in output_fasta:
        if i != 'NA':
            of.write(i)


