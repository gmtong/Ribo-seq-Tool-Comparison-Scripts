#This script takes the RiboTish_inputfasta.txt file and converts to a formatted fasta file
#Usage is python RiboTish_formattedfasta.py RiboTish_inputfasta.txt RiboTish_fastaformatted.fasta
import sys
import pandas as pd
import numpy as np

inputgtf  = sys.argv[1]
output_fasta = sys.argv[2]
outputfasta1 = open(output_fasta,'w')
RiboTish_gtf_formatted = pd.read_csv(inputgtf, sep="\"", header = None)

RiboTish_fasta = RiboTish_gtf_formatted[[5,3]]

RiboTish_fasta_file = ">" + RiboTish_fasta[5] + "\n" + RiboTish_fasta[3]
RiboTish_fasta_file = RiboTish_fasta_file.values.tolist()

RiboTish_fasta_cleaned = []
for i in RiboTish_fasta_file:
    removed = i.replace("*", "")
    RiboTish_fasta_cleaned.append(removed)

#output_file 
for k in RiboTish_fasta_cleaned:
   outputfasta1.write(str(k) + "\n")
outputfasta1.close()
