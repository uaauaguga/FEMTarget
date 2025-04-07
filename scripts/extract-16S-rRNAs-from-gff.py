#!/usr/bin/env python
import argparse
import sys
import os
import subprocess
from pyfaidx import Fasta
import logging
logging.basicConfig(level=logging.INFO, format='[%(asctime)s] [%(levelname)s] %(name)s: %(message)s')
logger = logging.getLogger('extract 16S rRNA')

def main():
    parser = argparse.ArgumentParser(description='Extract 16S rRNA')
    parser.add_argument('--gff',  '-g', type=str, required=True, help='Input genome annotation')
    parser.add_argument('--fasta','-f', type=str , help="Input genome sequence")
    parser.add_argument('--output','-o', type=str, required=True , help="Output rRNA sequences")
    args = parser.parse_args()
   

    logger.info("Extract rRNA coordinates ...")

    bed = args.gff + ".tmp.bed" 
    cmd = ["scripts/gff2bed.py", "--gff",  args.gff, "--bed", bed,  "--feature", "rRNA", "--name", "product"] 
    subprocess.run(cmd)

    logger.info("Extract rRNA sequences ...") 

    seq_id, start, end, strand = "", -1, -1, "."
    with open(bed) as f:
        for line in f:
            fields = line.strip().split("\t")
            if fields[3] == "16S ribosomal RNA":
                seq_id, start, end = fields[:3]
                strand = fields[5]
                start, end = int(start), int(end)
                break
    logger.info("Remove temp file ...")
    os.remove(bed)
    if start < 0:
        logger.warning("No annotated 16S rRNA")
    else:
        fout = open(args.output,"w")
        fasta = Fasta(args.fasta)
        sequence = fasta[seq_id][start:end]
        if strand == "-":
            sequence = sequence.reverse.complement
        sequence = str(sequence)
        print(f">{seq_id}:{start}-{end}({strand})",file=fout)
        print(sequence,file=fout)
        fout.close()
    logger.info("Finished.")               


if __name__ == "__main__":
    main()
