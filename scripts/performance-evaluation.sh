#!/bin/bash
asm_id=$1
gs=$2
params=$3
sf=$4
id=test/$asm_id/combined.final/${gs}.${params}
scripts/performance-evaluation.py -id $id  -kt benchmark/leader/$asm_id -o test.txt -sf "$sf" --min-homologs 0 -si tuning-test-sRNA-ids.txt #--validated
