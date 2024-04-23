#For all other tools besides RibORF0.1 for formatting fasta lists for tool comparison
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
            split2 = newline.split("chr")
            strand_split = split2[0]
            strand.append(strand_split[-1])
            split3 = split2[1]
            transcriptid.append(split3)
        else:
            newline1=i.replace('\n','')
            aasequence.append(newline1)
formattedids = pd.DataFrame({'strand':strand, 'transcriptid':transcriptid, 'aasequence':aasequence})
finalidformat = formattedids['strand'].astype(str) + formattedids['transcriptid'].astype(str) + '_' + formattedids['aasequence'].astype(str)
finalidformat = finalidformat.values.tolist()

for k in finalidformat:
   outfile1.write(str(k) + "\n")
outfile1.close()
