#Script for formatting the fasta files ouputted from RibORF0.1 fastas for tool comparison.
import sys
import numpy as np
import pandas as pd
inputfile  = sys.argv[1]
output_formattedgtf = sys.argv[2]
outfile1 = open(output_formattedgtf,'w')
transcriptid = []
aasequence = []
strand = []
with open(inputfile,"r") as f:
    lines = f.readlines()
    for i in lines:
        if i.startswith('>'):
            newline = i.replace('\n','')
            split = newline.split('_')
            id = split[0]
            split2 = id.split("chr")
            strand_split = split2[0]
            strand.append(strand_split[-1])
            split3 = split2[1]
            transcriptid.append(split3)
        else:
            newline1=i.replace('\n','')
            aasequence.append(newline1)
formattedids = pd.DataFrame({'strand':strand, 'transcriptid':transcriptid, 'aasequence':aasequence})
riborfidformat = formattedids['strand'].astype(str) + formattedids['transcriptid'].astype(str) + '_' + formattedids['aasequence'].astype(str)
riborfidformat = riborfidformat.values.tolist()

for k in riborfidformat:
   output_file.write(str(k) + "\n")
output_file.close()
