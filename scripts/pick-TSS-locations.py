#!/usr/bin/env python
import argparse
import subprocess
import os
import logging
import numpy as np
from collections import defaultdict
from scipy.stats import norm
import pandas as pd
from scipy.stats import gumbel_r
import io
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('pick TSS location')


def main():
    parser = argparse.ArgumentParser(description='prepare leader sequence location')
    parser.add_argument('--input','-i',type=str, required=True, help="score of promoter prediction in bed format")
    parser.add_argument('--gene','-g',type=str, required=True, help='gene interval in bed format, strandness is required')    
    parser.add_argument('--output','-o',type=str, required=True, help="where to save the scores")
    parser.add_argument('--contig', '-c', type=str, required=True, help="contig length")
    args = parser.parse_args()
    
    lengths = {}
    logger.info("Load contig length ...")
    with open(args.contig) as f:
        for line in f:
            seq_id, length = line.strip().split("\t")[:2]
            length = int(length)
            lengths[seq_id] = length
            
    
    logger.info("Load TSS scores ...")
    fwd_scores = {}
    rev_scores = {}
    with open(args.input) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, start, end = fields[:3]
            start, end = int(start), int(end)
            strand = fields[5]
            if seq_id not in fwd_scores:
                fwd_scores[seq_id] = np.zeros(lengths[seq_id])
                rev_scores[seq_id] = np.zeros(lengths[seq_id])
            score = float(fields[4])                
            if strand == "+":
                fwd_scores[seq_id][start] = score
            else:
                rev_scores[seq_id][start] = score
                
    logger.info("Impute offset length for each gene ...")   
    gene_id2iv = {}
    gene_id2leader_length = {}
    with open(args.gene) as f:
        for line in f:
            fields = line.strip().split("\t")
            seq_id, s, e = fields[:3]
            protein_id, strand = fields[3], fields[5]
            s, e = int(s), int(e)
            strand = fields[5] 
            gene_id2iv[protein_id] = (seq_id, s, e, strand)
            if s > e:
                print(*fields)
                continue
            if strand == "+":
                start = max(s-200,0)
                end = s
                if start > end:
                    start, end = end-1, start
                scores = fwd_scores[seq_id][start:end]
            else:
                start = e
                end =  min(lengths[seq_id],e+200)
                if start > end:
                    start, end = end-1, start
                scores = rev_scores[seq_id][start:end]
            positions = np.where(scores>0.7)[0]
            if len(positions) > 0:
                length = positions - 200
                gene_id2leader_length[protein_id] = length
    
    fout = open(args.output,"w")
    putative_length = sorted(list(gene_id2leader_length.values()))[int(0.9*length(gene_id2leader_length))]
    for protein_id in gene_id2iv:
        seq_id, s, e, strand = gene_id2iv[protein_id] 
        length = gene_id2leader_length.get(protein_id,putative_length)
        if strand == "+":
            start = max(s-length,0)
            end = min(lengths[seq_id], s + 100)
        else:
            start = max(e-100,0)
            end = min(lengths[seq_id], e + length)
         
        
        
    


if __name__ == "__main__":
    main()
    
            
