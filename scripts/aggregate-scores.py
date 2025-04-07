#!/usr/bin/env python
import logging
import argparse
import pandas as pd
import json
import re
from collections import defaultdict
import numpy as np
from scipy.stats import genextreme, norm, gumbel_r
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("score aggregation")

def get_weight(tree):
    tree = iter(tree)
    stack = []
    weights = defaultdict(float)
    eoc = False
    for c in tree:
        if c == ";":
            break
        if c  == " " or c == "\n":
            continue
        line = c
        while c not in "(),":
            c = next(tree)
            if c  == " " or c == "\n":
                continue
            line += c
        if line == "(":
            
            # start of a clade
            stack.append([])
        else:
            assert ":" in line
            #if line.startswith(":"):
            if eoc:
                d = float(line.split(":")[1][:-1])
                eoc = False
            else:
                name, d = line.split(":")
                d = float(d[:-1])
                name = [name]
            stack[-1].append((name,d))
            if line.endswith(")"):
                eoc = True
                # end of a clade
                attrs = stack.pop()
                # either leave or an internal node
                name = []
                for n ,d in attrs:
                    #print(n,d)
                    name.append(n)
                    leaves = n
                    for leave in leaves:
                        weights[leave] += d/len(leaves)
                flatened_name = []
                for name_list in name:
                    if isinstance(name_list,list):
                        flatened_name += name_list
                    else:
                        flatened_name += [name_list]
                name = flatened_name
    assert len(stack) == 0
    return weights

def get_combined_score(scores):
    xs = []
    ws = []
    for seq_id in scores:
        genome_id = seq_id.split(":")[0]
        xs.append(scores[seq_id])
        ws.append(genome_weights[genome_id])
    xs, ws = np.array(xs), np.array(ws)
    ws = ws/ws.sum()
    combined_score = (xs*ws).sum()/np.sqrt((1-correlation)*(ws*ws).sum() + correlation*ws.sum()**2)
    return combined_score

def main():
    parser = argparse.ArgumentParser(description='aggregate scores')
    parser.add_argument('--tree', '-t', type=str, required=True, help='tree file')
    parser.add_argument('--input', '-i', type=str, required=True, help='predicted  interaction')
    parser.add_argument('--correlation', '-c', type=float, default=0.8, help='correlation')
    parser.add_argument('--output', '-o', type=str, required=True, help='output statistics')
    parser.add_argument('--normalize','-n', action = "store_true" , help="whether normalize the score")
    parser.add_argument('--weights','-w', required=True, help="The weights for scoring")
    args = parser.parse_args()


    weights = json.load(open(args.weights))
    scores = defaultdict(dict)
    logger.info("Load scores ...")
    table = pd.read_csv(args.input,sep="\t")

    table["score"] = 0
    for key in weights:
        if weights[key] == 0:
            logger.info(f"{key} has zero weight, skip it.")
            continue
        assert key in table.columns, f"{key} not present in input table"
        values = table[key].fillna(0)
        # clip outliers
        values[values>7] = 7
        values[values<-7] = -7
        table["score"] += values*weights[key]
    if args.normalize:
        mean = table["score"].mean()
        std = table["score"].std()
        table["score"] = (table["score"] - mean)/std
        logger.info(f"Score mean is {mean}, stderr is {std}.")
    for query_id, target_id, score in table.loc[:,["query id", "protein id", "score"]].to_records(index=False):
        scores[query_id][target_id] = score

    global correlation, genome_weights
    correlation = args.correlation
    
    logger.info("extract sequence weight based on specified tree ...")
    global tree
    tree = open(args.tree).read()
    genome_weights = get_weight(tree)

    fout = open(args.output,"w")
    print("query id","target id","raw","denoised",sep="\t",file=fout)
    for query_id in scores:
        pair_scores = scores[query_id]
        combined_score = get_combined_score(pair_scores)
        for seq_id in pair_scores:
            print(query_id, seq_id, pair_scores[seq_id], combined_score, sep="\t", file=fout)
    fout.close()
    logger.info("All done.")
            
  
if __name__ == "__main__":
    main()
