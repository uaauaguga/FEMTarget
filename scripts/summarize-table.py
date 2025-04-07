#!/usr/bin/env python
import argparse
import os
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(name)s] [%(levelname)s] %(message)s')
logger = logging.getLogger("summarize tables")

def main():
    parser = argparse.ArgumentParser(description='Summarize several tables into a single matrix')
    parser.add_argument('--indir', '-i', type=str, required=True, help='directory contains input tables')
    parser.add_argument('--sample-ids', '-s', type=str, help='sample ids to include in the output matrix')
    parser.add_argument('--value-field', '-vf', type=int, required=True,help='value field in the table to summarize')
    parser.add_argument('--row-field', '-rf', type=int, required=True, help='row field in the table to summarize')
    parser.add_argument('--row-name', '-rn', type=str, default="feature", help='index name in the output matrix')
    parser.add_argument('--comment', '-c', default="#", type=str, help='line starts with this will be skipped')
    parser.add_argument('--first-line', '-f', type=int, default=0, help='only consider records after this line. all lines will be considered by default ')
    parser.add_argument('--output','-o', type=str, required=True, help='ouput matrix')
    parser.add_argument('--fillna', type=float, help='if specificied, fill na values with a float number')
    parser.add_argument('--integer', '-int', action = "store_true", default=False, help='whether output values should be interger')
    parser.add_argument('--sparse', '-sp', action = "store_true", default=False, help='whether write output as three column format (feature_id, sample_id, value). Useful for sparse feature.')
    parser.add_argument('--formatter',type=str,default="{}.txt", help='mapping from sample id to input file name')
    args = parser.parse_args()
    records = []
    if args.sample_ids is not None:
        sample_ids = open(args.sample_ids).read().strip().split("\n")
        print(sample_ids)
        logger.info("Check inputs ...")
        not_presented = []
        for sample_id in sample_ids:
            path = os.path.join(args.indir, args.formatter.format(sample_id))
            if not os.path.exists(path):
                not_presented.append(sample_id)
        if len(not_presented):
            logger.warning(f"{','.join(not_presented)} does not exists in input directory .")
        sample_ids = [x for x in sample_ids if x not in not_presented]
    else:
        sample_ids = []
        for txt in os.listdir(args.indir):
            if not txt.endswith(".txt"):
                continue
            sample_id = txt[:-4]
            sample_ids.append(sample_id)
        #sys.exit(1)
    logger.info("Load tables ...")
    for sample_id in tqdm(sample_ids):
        path = os.path.join(args.indir, args.formatter.format(sample_id))
        print(path)
        encountered_feature_ids = set()
        n_duplicated = 0
        with open(path) as f:
            i = 0
            if args.first_line > 0:
                for xx in range(args.first_line):
                    i += 1
                    _ = next(f)
            for line in f:
                i += 1
                line = line[:-1]
                if line.startswith(args.comment):
                    continue
                fields = line.split("\t")
                assert args.row_field < len(fields) and args.value_field < len(fields), f"in line {i}, only {len(fields)} available"                
                feature_id, value = fields[args.row_field], fields[args.value_field]
                if len(feature_id) == 0:
                    feature_id = "unassigned"
                if feature_id not in encountered_feature_ids:
                    if len(value) == 0:
                        value = np.nan
                    else:
                        value = round(float(value),3)
                    if value == 0:
                        continue
                    records.append((feature_id, sample_id, value))
                    encountered_feature_ids.add(feature_id)
                else:
                    print(feature_id)
                    n_duplicated += 1
            if n_duplicated > 0:
                logger.warning(f"Found {n_duplicated} duplicated features in {sample_id}, only consider the first one. ")
    table = pd.DataFrame.from_records(records)
    table.columns = [args.row_name,"sample_id","values"]
    if args.sparse:
        logger.info("Save three column matrix as sparse is specified ...")
        if args.integer:
            table["values"] = table["values"].astype(int)
        logger.info(f"Saving feature matrix to {args.output}")
        table.to_csv(args.output,sep="\t",index=False)
    else:
        logger.info("Pivot the table ...")
        matrix = table.pivot(index=args.row_name,columns="sample_id",values="values")
        if args.fillna is not None:
            matrix = matrix.fillna(args.fillna)
        if args.integer:
            matrix = matrix.fillna(0).astype(int)
        logger.info(f"Save feature matrix to {args.output} ...")
        matrix.to_csv(args.output,sep="\t") 
    logger.info("All done .")


if __name__ == "__main__":
    main()
