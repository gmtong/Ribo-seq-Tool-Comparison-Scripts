#From the ORFQuant .Rdata file, go into the R object and take the list of smORFS identified by ORFquant and create a csv file containing the transcript ID name, seqnames, start, end, length, and strand.
#once the CSV file has been created, this script takes the ORFquant_fasta output file and the csv file to create a fasta file.

#usage: input1:csv with names + genomic coordinates input2: fasta from orfquant output name. python3 input1 input2 output1
#python ORFquant_fasta_format.py ORFquant_csv.csv ORFquant_fastaoutput.fasta ORFquant_fasta_blastfiltered.fasta


import sys
import pandas as pd
import numpy as np
inputfile1 = sys.argv[1]

orfquant_input = pd.read_csv(inputfile1, sep=",", header='infer')
orfquant_input['seq_strand_seqnames'] = orfquant_input['x'].astype(str) + "_" + orfquant_input['strand'].astype(str) + "_" + orfquant_input['seqnames'].astype(str)
orfquant_input['start_end'] = orfquant_input['start'].astype(str) + " " + orfquant_input['end'].astype(str)
orfquant_input.drop("width", axis=1,inplace=True)
orfquant_input.drop("strand", axis=1,inplace=True)
orfquant_input.drop("x", axis=1,inplace=True)
orfquant_input.drop("seqnames", axis=1,inplace=True)
orfquant_input.drop("start", axis=1,inplace=True)
orfquant_input.drop("end", axis=1,inplace=True)


orfquant_input = orfquant_input.groupby(['seq_strand_seqnames']).start_end.apply(' '.join).reset_index(name='listofcoordinates')


listofcoordinates = orfquant_input['listofcoordinates']
listofcoordinates = listofcoordinates.values.tolist()
ids = orfquant_input['seq_strand_seqnames']
ids = ids.values.tolist()
ids_df = pd.DataFrame({'ids':ids})
ids_df_split = ids_df['ids'].str.split('\_', expand = True)


start = []
stop = []

for i in listofcoordinates:
    strip1 = i.strip().split(' ')
    start.append(strip1[0])
    stop.append(strip1[-1])

genomiccoordinates_df = pd.DataFrame({'start': start, 'stop':stop})


df_format = pd.concat([ids_df_split, genomiccoordinates_df], axis=1, join='inner')

transcript_relativeid = df_format[0] + "_" +  df_format[1] + "_" + df_format[2]
transcript_genomicid =  df_format[0] + df_format[3] + df_format[4] + ":" + df_format['start'] + "-" + df_format['stop']

genomic_coord_df_match = pd.DataFrame({'transcript_relative': transcript_relativeid, 'transcript_genomic': transcript_genomicid})

add_delimiter = []
inputfile2 = sys.argv[2]

with open(inputfile2, 'r') as f:
    lines = f.readlines()
    for i in lines:
        if i.startswith('>'):
            replace = i.replace('\n', '|')
            add_delimiter.append(replace)
        else:
            add_delimiter.append(i)
newlines_removed = []
for i in add_delimiter:
    strip1 = i.replace('\n','')
    newlines_removed.append(strip1)
joined_format = ''.join(newlines_removed)
joined_format = joined_format.split('>')

            
fasta_df = pd.DataFrame({'fasta': joined_format})
fasta_df.drop(index=fasta_df.index[0], axis = 0, inplace=True)
fasta_split_match = fasta_df['fasta'].str.split("\|", expand = True)
fasta_split_match.columns = ['transcript_relative', 'annotation_type', 'transcript_id', 'starT_stop_type', 'annotation', 'AAsequence']
#print(fasta_split_match)

merged_df = pd.merge(genomic_coord_df_match, fasta_split_match, on = 'transcript_relative')
filteredmerged_df = merged_df[merged_df["starT_stop_type"].str.contains("novel_Internal")==False]
#print(len(filteredmerged_df))
#print(len(merged_df))

extract_genomic_df = '>' + filteredmerged_df['transcript_genomic'] + '\n' + filteredmerged_df['AAsequence']
extract_genomic_df = extract_genomic_df.values.tolist()

out = sys.argv[3]

output_file = open(out, 'w')
for k in extract_genomic_df:
   output_file.write(str(k) + "\n")
output_file.close()

#After running this script, the outputted fasta was run in a blastp filter following our inhouse pipeline to look for unannotated smORFs.
#the outputted fasta will be used in the next script.
