#This script is to format the RiboCode output.txt file to create a fasta file identified by RiboCode
#Usage: python RiboCode_fasta_formatting.py Input_RiboCode_output.txt RiboCode_fasta_formatted.fa

import sys
import pandas as pd
inputgtf  = sys.argv[1]
output_formattedgtf = sys.argv[2]
outfile1 = open(output_formattedgtf,'w')
Ribocode_gtf = pd.read_csv(inputgtf, sep="\t", header='infer')
Ribocode_Extraction = Ribocode_gtf[['ORF_ID', 'transcript_id', 'strand', 'chrom','AAseq']]

positivestrand_ribocode_extraction = Ribocode_Extraction[Ribocode_Extraction.strand == '+']
negativestrand_ribocode_extraction = Ribocode_Extraction[Ribocode_Extraction.strand == '-']

pos_strand_ribocode_extraction = positivestrand_ribocode_extraction['ORF_ID'].str.split("\_", expand = True)
pos_strand_ribocode_extraction.columns = ['Ensembl_ID', 'Start', 'Stop', 'Extra']
pos_strand_ribocode_extraction['Stop'] = pos_strand_ribocode_extraction['Stop'].astype(int) - 3 

Ribocode_pos_fasta = ">" + positivestrand_ribocode_extraction['transcript_id'] + positivestrand_ribocode_extraction['strand'] + positivestrand_ribocode_extraction['chrom'] + ":" + pos_strand_ribocode_extraction['Start'].astype(str) + "-" + pos_strand_ribocode_extraction['Stop'].astype(str) + "\n" + positivestrand_ribocode_extraction['AAseq']


neg_strand_ribocode_extraction = negativestrand_ribocode_extraction['ORF_ID'].str.split("\_", expand = True)
neg_strand_ribocode_extraction.columns = ['Ensembl_ID', 'Start', 'Stop', 'Extra']
neg_strand_ribocode_extraction['Stop'] = neg_strand_ribocode_extraction['Stop'].astype(int) + 3 

Ribocode_neg_fasta = ">" + negativestrand_ribocode_extraction['transcript_id'] + negativestrand_ribocode_extraction['strand'] + negativestrand_ribocode_extraction['chrom'] + ":" + neg_strand_ribocode_extraction['Stop'].astype(str) + "-" + neg_strand_ribocode_extraction['Start'].astype(str) + "\n" + negativestrand_ribocode_extraction['AAseq']


Ribocode_pos_fasta = Ribocode_pos_fasta.values.tolist()
Ribocode_neg_fasta = Ribocode_neg_fasta.values.tolist()

Ribocode_fasta_joined = Ribocode_pos_fasta + Ribocode_neg_fasta

for k in Ribocode_fasta_joined:
   outfile1.write(str(k) + "\n")
outfile1.close()
