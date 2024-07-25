#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 10:04:52 2023

@author: nathan
"""


"""
Creating BED files that work with the pseudoref the way I made it.

Basically, every ID is treated as its own "chromosome" with an independent index so it casused things like bedtools to not play nice.

I am writing this python script to generate BED files for bedtools coverage for the TE junctions.
"""
import argparse

def parse_args():
    parser = argparse.ArgumentParser(\
                    description='Creating a bed file for 3bp encompassing the left and right junctions of the SV and the 3bp encompassing the polymorphic site/insertion point.')
    parser.add_argument('-f', '--fasta_index', action='store',\
                        type=str,required=True,help='Pseudoreference fasta index generated from "samtools faidx". Should have IDs for all SV-present alleles and SV-absent alleles/insertion points')
    parser.add_argument('-o', '--out_file', action='store',\
                        type=str,help='Write output to outfile')
                        
    args = parser.parse_args()
    return args

def insertion_bed(line3):
    start = 301-2 # insertion point is at 1-based position 301
    end = 301 + 1
    ins_l = [line3[0],str(start),str(end),'insertion_point']
    return ins_l

    
def left_right_sv_juncs_bed(line3):
    lstart = 300-2
    lend = 300+1
    rstart = int(float((line3[1])) - float(300) - float(2))
    rend = int(float((line3[1])) - float(300) +float(1))
    left_junc_l = [line3[0],str(lstart),str(lend), 'left']
    right_junc_l = [line3[0],str(rstart),str(rend), 'right']
    return left_junc_l , right_junc_l

def write_lists_to_output(lists, output):
    for lst in lists:
        line = '\t'.join(map(str, lst))
        output.write(line + '\n')


def main():
    
    args=parse_args()
    file = args.fasta_index
    lines_to_print = list()
    with open(file, 'r') as fh:
        for line in fh:
            line2 = line.strip()
            line3 = line2.split()
            if "ins" in line3[0]:
                left_junc_l, right_junc_l = left_right_sv_juncs_bed(line3)
                lines_to_print.append(left_junc_l)
                lines_to_print.append(right_junc_l)
            else:
                ins_l = insertion_bed(line3)
                lines_to_print.append(ins_l)
    if args.out_file:
        with open(args.out_file, 'w') as of:
            write_lists_to_output(lines_to_print, of)
    else:
        import sys
        write_lists_to_output(lines_to_print, sys.stdout)
    

if __name__ == "__main__":
    main()