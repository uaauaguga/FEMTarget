#!/usr/bin/env python
import subprocess
import logging
import pandas as pd
import json
import argparse
import numpy as np
import io
import os
import re
from collections import defaultdict
import pandas as pd
from ete3 import Tree
import subprocess
from multiprocessing import Pool
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('comparative scoring')

def prune(tree, names):
    nodes = []
    for node in tree.traverse():
        name = node.name
        if name in names:
            nodes.append(node)
    tree.prune(nodes)
    return tree


def make_nex(path, scores):
    scores2 = defaultdict(dict)
    genome_ids = set()
    query_ids = sorted(list(scores.keys()))
    for query_id in scores:
        for target_id in scores[query_id]:
            genome_id = target_id.split(":")[0]
            scores2[genome_id][query_id] = scores[query_id][target_id]
            genome_ids.add(genome_id)
    scores = scores2
    genome_ids = sorted(list(genome_ids))
    ntax = len(genome_ids)
    nchar = len(query_ids)
    fout = open(path,"w")
    print("#NEXUS", file=fout)
    print("BEGIN DATA;",file=fout)
    print(f"DIMENSIONS  NCHAR={nchar} NTAX={ntax};", file=fout)
    print("FORMAT datatype=Continuous missing=? gap=-;",file=fout)
    print("",file=fout)
    print("    MATRIX",file=fout)
            
    for genome_id in genome_ids:
        scores_by_genome = [str(round(scores[genome_id][q],4)) for q in query_ids]
        scores_by_genome = " ".join(scores_by_genome)
        print(f"    {genome_id} {scores_by_genome}", file=fout)
    print("",file=fout)
    print(";",file=fout)
    print("END;",file=fout)
    fout.close()

def main():
    parser = argparse.ArgumentParser(description='scoring sRNA targt interaction with comparative scoring')
    parser.add_argument('--input',  '-i', type=str, required=True, help='interaction scores')
    parser.add_argument('--output','-o', type=str , required=True, help="where to save output")
    parser.add_argument('--select','-s', type=str , required = True, help="target ids to consider")
    parser.add_argument('--full-tree','-ft', type=str , required=True, help="the phylogenetic tree")
    parser.add_argument('--pruned-tree','-pt', type=str , required=True, help="the phylogenetic tree")
    parser.add_argument('--normalize','-n', action = "store_true" , help="whether normalize the score")
    parser.add_argument('--weights','-w', default="models/scoring.json", help="The weights for scoring")
    args = parser.parse_args()


    weights = json.load(open(args.weights))
    scores = defaultdict(dict)
    logger.info("Load scores ...")
    table = pd.read_csv(args.input,sep="\t")
    
    table["score"] = 0
    for key in weights:
        if key not in table.columns:
            logger.info(f"{key} not in input table. Skip it.") 
            continue
        values = table[key].fillna(0)
        values[values>4] = 4
        values[values<-4] = -4       
        table["score"] += values*weights[key]

    if args.normalize:
        mean = table["score"].mean()
        std = table["score"].std()
        table["score"] = (table["score"] - mean)/std    
        logger.info(f"Score mean is {mean}, stderr is {std}.")

    protein_ids = []
    logger.info("Load selected protein ids ...")
    with open(args.select) as f:
        for line in f:
            fields = line.strip().split("\t")
            protein_id = fields[0]
            protein_ids.append(protein_id)

    
    genome_ids = set()
    for query_id, target_id, score in table.loc[:,["query id", "protein id", "score"]].to_records(index=False):
        protein_id = query_id.split(":")[1].split("-")[0]
        genome_ids.add(target_id.split(":")[0])
        if protein_id not in protein_ids:
            continue
        scores[query_id][target_id] = score    
    n_targets = len(scores)
    n_genomes = len(genome_ids)
    logger.info(f"Select {n_targets} proteins for consideration.")
    logger.info(f"{n_genomes} genomes contains the sRNA")
    
    scores2 = defaultdict(dict)
    for query_id in scores:
        if len(scores[query_id]) == n_genomes:                              
            for target_id in scores[query_id]:
                scores2[query_id][target_id] = scores[query_id][target_id]
    scores = scores2
    n_targets = len(scores)
    logger.info(f"{n_targets} proteins are finally used.")

    logger.info("Load tree ...")
    full_tree = Tree(newick=args.full_tree,format=1,quoted_node_names=True)
    genome_ids = sorted(list(genome_ids))
    used_tree = prune(full_tree,genome_ids)
    nex = args.output

    logger.info("Make nex ...")
    make_nex(nex, scores)

    logger.info("Save pruned tree ...")
    with open(args.pruned_tree,"w") as f:
        f.write(used_tree.write(format=1))

    logger.info("All done.")


                

if __name__ == "__main__":
    main()
