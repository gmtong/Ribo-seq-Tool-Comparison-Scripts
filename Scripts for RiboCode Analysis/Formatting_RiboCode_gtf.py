#Outputs formatted gtf from RiboCode output, option to run sed 's/_PAR_Y//g' command to remove pseudoautosomal zones from the output files before running script
#Usage: python RiboCode_gtf_format.py Input_RiboCode_gtf Output_RiboCode_gtf_formatted.gtf

import sys 
import pandas as pd
import numpy as np

inputgtf  = sys.argv[1]
ribocode_gtf = pd.read_csv(inputgtf, sep="\t", header=None)
output_formattedgtf = sys.argv[2]
outfile1 = open(output_formattedgtf,'w')
ribocode_gtf.columns = ['chr', 'Tool', 'exon', 'Start', 'Stop', 'Dot1', 'Strand', 'Dot2', 'Transcript_info']

positivestrand_ribocode = ribocode_gtf[ribocode_gtf.Strand == '+']
negativestrand_ribocode = ribocode_gtf[ribocode_gtf.Strand == '-']

ribocode_gtf_positive_transcriptid_split = positivestrand_ribocode['Transcript_info'].str.split("\"", expand = True)

ribocode_positivestrand_transcriptid = ribocode_gtf_positive_transcriptid_split[1].str.split("_", expand = True)
ribocode_positivestrand_transcriptid.columns = ['geneid','Start', 'Stop', 'Extra']



ribocode_positivestrand_transcriptid['Stop'] = ribocode_positivestrand_transcriptid['Stop'].astype(int) - 3


ribocode_gtf_negative_transcriptid_split = negativestrand_ribocode['Transcript_info'].str.split("\"", expand = True)

ribocode_negativestrand_transcriptid = ribocode_gtf_negative_transcriptid_split[1].str.split("\_", expand = True)
ribocode_negativestrand_transcriptid.columns = ['geneid','Start', 'Stop', 'Extra']

ribocode_negativestrand_transcriptid['Stop'] = ribocode_negativestrand_transcriptid['Stop'].astype(int) + 3



positivestrand_corrected = positivestrand_ribocode['chr'] + "\t" + positivestrand_ribocode['Tool'] + "\t" + positivestrand_ribocode['exon'] + "\t" + positivestrand_ribocode['Start'].astype(str) + "\t" + positivestrand_ribocode['Stop'].astype(str) + "\t" + positivestrand_ribocode['Dot1'] + "\t" + positivestrand_ribocode['Strand'] + "\t" + positivestrand_ribocode['Dot2'] + "\t" + ribocode_gtf_positive_transcriptid_split[0] + " " + "\"" + ribocode_gtf_positive_transcriptid_split[1] + "\"" + ribocode_gtf_positive_transcriptid_split[2] + "\"" + ribocode_gtf_positive_transcriptid_split[3] + "\"" + ribocode_gtf_positive_transcriptid_split[4] + "\"" + ribocode_gtf_positive_transcriptid_split[9] + positivestrand_ribocode['Strand'] + positivestrand_ribocode['chr'] + ":" + ribocode_positivestrand_transcriptid['Start'].astype(str) + "-" + ribocode_positivestrand_transcriptid['Stop'].astype(str) + "\"" + ribocode_gtf_positive_transcriptid_split[6] + "\"" + ribocode_gtf_positive_transcriptid_split[7] + "\"" + ribocode_gtf_positive_transcriptid_split[8] + "\"" + ribocode_gtf_positive_transcriptid_split[9] + positivestrand_ribocode['Strand'] + positivestrand_ribocode['chr'] + ":" + ribocode_positivestrand_transcriptid['Start'].astype(str) + "-" + ribocode_positivestrand_transcriptid['Stop'].astype(str) + "\""
positivestrand_corrected = positivestrand_corrected.values.tolist()


negativestrand_corrected = negativestrand_ribocode['chr'] + "\t" + negativestrand_ribocode['Tool'] + "\t" + negativestrand_ribocode['exon'] + "\t" + negativestrand_ribocode['Start'].astype(str) + "\t" + negativestrand_ribocode['Stop'].astype(str) + "\t" + negativestrand_ribocode['Dot1'] + "\t" + negativestrand_ribocode['Strand'] + "\t" + negativestrand_ribocode['Dot2'] + "\t" + ribocode_gtf_negative_transcriptid_split[0] + " " + "\"" + ribocode_gtf_negative_transcriptid_split[1] + "\"" + ribocode_gtf_negative_transcriptid_split[2] + "\"" + ribocode_gtf_negative_transcriptid_split[3] + "\"" + ribocode_gtf_negative_transcriptid_split[4] + "\"" + ribocode_gtf_negative_transcriptid_split[9] + negativestrand_ribocode['Strand'] + negativestrand_ribocode['chr'] + ":" + ribocode_negativestrand_transcriptid['Stop'].astype(str) + "-" + ribocode_negativestrand_transcriptid['Start'].astype(str) + "\"" + ribocode_gtf_negative_transcriptid_split[6] + "\"" + ribocode_gtf_negative_transcriptid_split[7] + "\"" + ribocode_gtf_negative_transcriptid_split[8] + "\"" + ribocode_gtf_negative_transcriptid_split[9] + negativestrand_ribocode['Strand'] + negativestrand_ribocode['chr'] + ":" + ribocode_negativestrand_transcriptid['Stop'].astype(str) + "-" + ribocode_negativestrand_transcriptid['Start'].astype(str) + "\""
negativestrand_corrected = negativestrand_corrected.values.tolist()


joined_ribocode = positivestrand_corrected + negativestrand_corrected

for k in joined_ribocode:
   outfile1.write(str(k) + "\n")
outfile1.close()
