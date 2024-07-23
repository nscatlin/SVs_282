#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  5 14:03:42 2022

@author: nathan
"""

import pysam
import re
from collections import defaultdict
import argparse

"""
Usage:
    
    python anchorwave_include_sv_insertion_fastagen.py \
        anchorwave_tsv_file.tsv 
        ref_genotype \
        query_genotype \
        ref_genome.fasta \
        query_genome.fasta


"""
def parse_args():
    parser = argparse.ArgumentParser(\
                    description='Extract polymorphic regions surrounding SVs in SV-present and SV-absent genontypes.')
    parser.add_argument('-t', '--tsv_file', action='store',\
                        type=str,required=True,help='tab seperated file (TSV) denoting locations of polymorphic SVs between your ascertainment set')
    parser.add_argument('-a', '--genotype_1', action='store',\
                        type=str,required=True,help='One of your two genotypes present in your ascertainment set')
    parser.add_argument('-b', '--genotype_2', action='store',\
                        type=str,required=True,help='The other of the two genotypes present in your ascertainment set')
    parser.add_argument('-g', '--genome_1', action='store',\
                        type=str,required=True,help='Reference genome of "-a/--genotype_1" genoytpe (FASTA)')
    parser.add_argument('-j', '--genome_2', action='store',\
                        type=str,required=True,help='Reference genome of "-b/--genotype_2" genoytpe (FASTA)')
    args = parser.parse_args()
    return args


def tsv_processing(chr_tsv_file,genotype1, genotype2, genome1, genome2):
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
            if Block_Type.split('_in')[-1] == genotype1:
    
                # chromosome_num_ref = line3[0]
                ins_genotype = genotype1
                non_ins_genotype = genotype2
                ins_fasta = genome1
                non_ins_fasta = genome2
    
                
                ins_chromosome_num = line3[0]
                rowID = line3[12]
                # Block_Type = current_line_list[6]
                ins_start = line3[1]
                ins_end = line3[2]
                left_aln_region = line3[14]
                right_aln_region = line3[17]
                non_ins_coord = line3[4]
                non_ins_chromosome = line3[3]
    
                
                ins, non_ins = extract_sv_plus_adjacent_seqs(ins_genotype, 
                                                             non_ins_genotype,
                                                             ins_fasta, non_ins_fasta, 
                                                             ins_chromosome_num, rowID,
                                                             ins_start, ins_end,
                                                             non_ins_coord,
                                                             non_ins_chromosome, 
                                                             left_aln_region,
                                                             right_aln_region)
                output_fasta.extend([ins, non_ins])
    
            else:
    
                ins_genotype = genotype2
                non_ins_genotype = genotype1
                ins_fasta = genome2
                non_ins_fasta = genome1            
                ins_chromosome_num = line3[3]
                rowID = line3[12]
                # Block_Type = current_line_list[6]
                ins_start = line3[4]
                ins_end = line3[5]
                left_aln_region = line3[14]
                right_aln_region = line3[17]
                non_ins_coord = line3[1]
                non_ins_chromosome = line3[0]            
                ins, non_ins = extract_sv_plus_adjacent_seqs(ins_genotype,
                                                             non_ins_genotype,
                                                             ins_fasta,
                                                             non_ins_fasta,
                                                             ins_chromosome_num,
                                                             rowID, 
                                                             ins_start,
                                                             ins_end, 
                                                             non_ins_coord, 
                                                             non_ins_chromosome,
                                                             left_aln_region,
                                                             right_aln_region)
    
                output_fasta.extend([ins, non_ins])
    return output_fasta

def extract_sv_plus_adjacent_seqs(ins_geno, non_ins_genotype, 
                                  insertion_fasta_file, absent_insertion_fasta_file, 
                                  ins_chromosome_num, rowID, ins_start, ins_end,
                                  non_ins_coord, non_ins_chromosome, left_aln_region, right_aln_region):
    """
    Extracts SV + 300bp flanking alignable regions in insertion genotype and 300 bp flanking alignable region in opposite genotype where insertion is absent
    """
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

        """
        Extact 300 bp flanking alignable region at insertion location in genome WITHOUT insertion
        """

        non_ins_start = int(float(non_ins_coord) - 300)
        non_ins_end = int(float(non_ins_coord) + 300)
        non_ins_fasta = pysam.faidx(absent_insertion_fasta_file, non_ins_chromosome+':'+str(non_ins_start)+'-'+str(non_ins_end))
        non_ins_fasta_re_change_fasta_id = re.sub(r'\>\S+', '>'+non_ins_genotype+'_chr'+non_ins_chromosome+':'+str(
            non_ins_start)+'-'+str(non_ins_end)+'_'+rowID, non_ins_fasta)
        return ins_fasta_re_change_fasta_id, non_ins_fasta_re_change_fasta_id
    else:
        ins_fasta_re_change_fasta_id = 'NA'
        non_ins_fasta_re_change_fasta_id = 'NA'
        return ins_fasta_re_change_fasta_id, non_ins_fasta_re_change_fasta_id




#################################################################################################################################################################################################################
#################################################################################################################################################################################################################
#################################################################################################################################################################################################################


def main():
    args=parse_args()
    chr_tsv_file = args.tsv_file
    genotype1 = args.genotype_1
    genotype2 = args.genotype_2
    genome1 = args.genome_1
    genome2 = args.genome_2
    
    output_fasta = tsv_processing(chr_tsv_file,genotype1, genotype2, genome1, genome2)
    

    
    with open(genotype1 + '_' +genotype2 +'_pseudoref.fasta', 'w') as of:
        for i in output_fasta:
            if i != 'NA':
                of.write(i)

if __name__ == "__main__":
    main()




