#module load ucsc-tools and use genePredToGtf to create a gtf file. Then run the following script
#usage: python python_genepredtogtf.gtf genepredtogtf_output.gtf gtf_formatted.gtf
#following script will format the gtf file.

import sys
import numpy as np
import pandas as pd
pd.options.mode.chained_assignment = None
inputgtf = sys.argv[1]
out = sys.argv[2]

#input file is the gtf file from genepredtogtf.gtf
genepred_gtf = pd.read_csv(inputgtf, sep='\t', header = None)
genepred_gtf.columns = ['chr', 'RibORF1.0', 'type', 'start', 'stop', 'dot1', 'strand', 'dot2', 'gene_id']
genepred_gtf['type'] = genepred_gtf['type'].str.strip()
genepred_gtf = genepred_gtf[genepred_gtf['type'].str.contains("exon")==False]
genepred_gtf = genepred_gtf[genepred_gtf['type'].str.contains("start_codon")==False]
genepred_gtf2 = genepred_gtf[genepred_gtf['type'].str.contains("stop_codon")==False]
genepred_gtf2 = genepred_gtf2[genepred_gtf2['type'].str.contains("transcript")==False]
genepred_gtf2['type'] = genepred_gtf2['type'].replace('CDS', 'exon')

genepred_format = genepred_gtf2['chr'].astype(str) + "\t" + genepred_gtf2['RibORF1.0'].astype(str) + "\t" + genepred_gtf2['type'].astype(str) + "\t" + genepred_gtf2['start'].astype(str) + "\t" + genepred_gtf2['stop'].astype(str) + "\t" + genepred_gtf2['dot1'].astype(str) + "\t" + genepred_gtf2['strand'].astype(str) + "\t" + genepred_gtf2['dot2'].astype(str) + "\t" + genepred_gtf2['gene_id'].astype(str)
genepred_format = genepred_format.values.tolist()

#output gtf formatted
output_file = open(out, 'w')
for j in genepred_format :
    output_file.write(str(j) + "\n")
output_file.close
