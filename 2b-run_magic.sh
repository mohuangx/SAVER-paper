#!/bin/bash
# Run MAGIC on datasets

# Melanoma Drop-seq

MAGIC.py -d SAVER-data/melanoma_dropseq.csv -o SAVER-data/melanoma_dropseq_magic.csv --cell-axis 'columns' csv 

# Baron 
MAGIC.py -d SAVER-data/baron_human_samp.csv -o SAVER-data/baron_human_samp_magic.csv --cell-axis 'columns' csv 

# Chen
MAGIC.py -d SAVER-data/chen_samp.csv -o SAVER-data/chen_samp_magic.csv --cell-axis 'columns' csv 

# La Manno
MAGIC.py -d SAVER-data/manno_human_samp.csv -o SAVER-data/manno_human_samp_magic.csv --cell-axis 'columns' csv 

# Zeisel
MAGIC.py -d SAVER-data/zeisel_samp.csv -o SAVER-data/zeisel_samp_magic.csv --cell-axis 'columns' csv 

# Clean results
Rscript bin/clean_magic.R

rm *magic.csv
