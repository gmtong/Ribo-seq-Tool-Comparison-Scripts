#RiboTish takes input RiboTishoutput.txt and outputs a fastainputfile and genepred text file
#Usage: python RiboTish_input_convert_Fasta_genepred.py RiboTishoutput.txt RiboTish_inputfasta.txt RiboTishgenepred.genepred
import sys
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
inputgtf  = sys.argv[1]

RiboTish_gtf = pd.read_csv(inputgtf, sep="\t", header='infer')
RiboTish_gtf = RiboTish_gtf.loc[RiboTish_gtf['AASeq'].str.len() < 151]

RiboTish_GenomePos = RiboTish_gtf['GenomePos'].str.split(":", expand = True)
RiboTish_GenomePos.columns = ['Chr', 'Start_Stop', 'Strand']

RiboTish_Start_Stop = RiboTish_GenomePos['Start_Stop'].str.split("-", expand = True)
RiboTish_Start_Stop.columns = ['Start', 'Stop']

RiboTish_gtf_new = RiboTish_GenomePos['Chr'] + "\t" + "RiboTish" + "\t" + "ORF" + "\t" + RiboTish_Start_Stop['Start'] + "\t" + RiboTish_Start_Stop['Stop'] + "\t" + "." + "\t" + RiboTish_GenomePos['Strand'] + "\t" + "." + "\t" + "transcript_id " + '"' + RiboTish_gtf["Tid"] + '"; ' + '"' + RiboTish_gtf['AASeq'] + '"' + "\t" + RiboTish_gtf["Blocks"]
RiboTish_gtf_new = RiboTish_gtf_new.str.split("\t", expand = True)
RiboTish_gtf_new.columns = ['chr', 'Tool', 'ORF', 'Start', 'Stop', 'Dot1', 'Strand', 'Dot2', 'TranscriptInfo', 'Blocks']


positivestrand_ribotish = RiboTish_gtf_new[RiboTish_gtf_new.Strand == '+']
negativestrand_ribotish = RiboTish_gtf_new[RiboTish_gtf_new.Strand == '-']

positivestrand_ribotish['Start'] = positivestrand_ribotish['Start'].astype(int) + 1
positivestrand_ribotish['Stop'] = positivestrand_ribotish['Stop'].astype(int) - 3
positivestrand_ribotish_transcriptinfo = positivestrand_ribotish['TranscriptInfo'].str.split("\"", expand = True)
positivestrand_ribotish_formatted = positivestrand_ribotish['chr'] + "\t" + positivestrand_ribotish['Tool'] + "\t" + 'transcript' + "\t"+ positivestrand_ribotish['Start'].astype(str) + "\t" + positivestrand_ribotish['Stop'].astype(str) + "\t" + positivestrand_ribotish['Dot1'] + "\t" + positivestrand_ribotish['Strand'] + "\t"+ positivestrand_ribotish['Dot2'] + "\t" + "gene_id" + " " + "\"" + positivestrand_ribotish_transcriptinfo[1] + positivestrand_ribotish['Strand'] + positivestrand_ribotish['chr'] + ":" + positivestrand_ribotish['Start'].astype(str) + "-" + positivestrand_ribotish['Stop'].astype(str) + "\"; "  + "AASeq" + " \"" + positivestrand_ribotish_transcriptinfo[3] + "\"; " + positivestrand_ribotish_transcriptinfo[0] + "\"" + positivestrand_ribotish_transcriptinfo[1] + positivestrand_ribotish['Strand'] + positivestrand_ribotish['chr'] + ":" + positivestrand_ribotish['Start'].astype(str) + "-" + positivestrand_ribotish['Stop'].astype(str) + "\"" + ";" + "\t" + positivestrand_ribotish['Blocks']


negativestrand_ribotish['Start'] = negativestrand_ribotish['Start'].astype(int) + 4
negativestrand_ribotish['Stop'] = negativestrand_ribotish['Stop'].astype(int)
negativestrand_ribotish_transcriptinfo = negativestrand_ribotish['TranscriptInfo'].str.split("\"", expand = True)
negativestrand_ribotish_formatted = negativestrand_ribotish['chr'] + "\t" + negativestrand_ribotish['Tool'] + "\t" + 'transcript' + "\t"+ negativestrand_ribotish['Start'].astype(str) + "\t" + negativestrand_ribotish['Stop'].astype(str) + "\t" + negativestrand_ribotish['Dot1'] + "\t" + negativestrand_ribotish['Strand'] + "\t"+ negativestrand_ribotish['Dot2'] + "\t" + "gene_id" + " " + "\"" + negativestrand_ribotish_transcriptinfo[1] + negativestrand_ribotish['Strand'] + negativestrand_ribotish['chr'] + ":" + negativestrand_ribotish['Start'].astype(str) + "-" + negativestrand_ribotish['Stop'].astype(str) + "\"; "  + "AAseq" + " \"" + negativestrand_ribotish_transcriptinfo[3] + "\"; " + negativestrand_ribotish_transcriptinfo[0] + "\"" + negativestrand_ribotish_transcriptinfo[1] + negativestrand_ribotish['Strand'] + negativestrand_ribotish['chr'] + ":" + negativestrand_ribotish['Start'].astype(str) + "-" + negativestrand_ribotish['Stop'].astype(str) + "\"" + ";" + "\t" + negativestrand_ribotish['Blocks']


positivestrand_ribotish_formatted = positivestrand_ribotish_formatted.values.tolist()
negativestrand_ribotish_formatted = negativestrand_ribotish_formatted.values.tolist()

joined_ribotish_gtf = positivestrand_ribotish_formatted + negativestrand_ribotish_formatted

joinedlist = pd.DataFrame(joined_ribotish_gtf).astype(str)
joinedlist = joinedlist[0].str.split("\t", expand = True)


#taking the exon blocks and splitting them (might need to adjust in the future if you have exons spanning longer than 10 exons)
joinedlist_blocks = joinedlist[9].str.split(",", expand = True)
print(joinedlist_blocks)
joinedlist_exon1 = joinedlist_blocks[0].str.split("-", expand = True)
joinedlist_exon1.columns = ['start_exon1', 'stop_exon1']
joinedlist_exon2 = joinedlist_blocks[1].str.split("-", expand = True)
joinedlist_exon2.columns = ['start_exon2', 'stop_exon2']
joinedlist_exon3 = joinedlist_blocks[2].str.split("-", expand = True)
joinedlist_exon3.columns = ['start_exon3', 'stop_exon3']
joinedlist_exon4 = joinedlist_blocks[3].str.split("-", expand = True)
joinedlist_exon4.columns = ['start_exon4', 'stop_exon4']
joinedlist_exon5 = joinedlist_blocks[4].str.split("-", expand = True)
joinedlist_exon5.columns = ['start_exon5', 'stop_exon5']
joinedlist_exon6 = joinedlist_blocks[5].str.split("-", expand = True)
joinedlist_exon6.columns = ['start_exon6', 'stop_exon6']
joinedlist_exon7 = joinedlist_blocks[6].str.split("-", expand = True)
joinedlist_exon7.columns = ['start_exon7', 'stop_exon7']
joinedlist_exon8 = joinedlist_blocks[7].str.split("-", expand = True)
joinedlist_exon8.columns = ['start_exon8', 'stop_exon8']
joinedlist_exon9 = joinedlist_blocks[8].str.split("-", expand = True)
joinedlist_exon9.columns = ['start_exon9', 'stop_exon9']
joinedlist_exon10 = joinedlist_blocks[9].str.split("-", expand = True)
joinedlist_exon10.columns = ['start_exon10', 'stop_exon10']

exon_starts = joinedlist_exon1['start_exon1'].astype(str) + "," + joinedlist_exon2['start_exon2'].astype(str) + "," + joinedlist_exon3['start_exon3'].astype(str) + "," + joinedlist_exon4['start_exon4'].astype(str) + "," + joinedlist_exon5['start_exon5'].astype(str) + "," + joinedlist_exon6['start_exon6'].astype(str) + "," + joinedlist_exon7['start_exon7'].astype(str) + "," + joinedlist_exon8['start_exon8'].astype(str) + "," + joinedlist_exon9['start_exon9'].astype(str) + "," + joinedlist_exon10['start_exon10'].astype(str) + ","
exon_stops =  joinedlist_exon1['stop_exon1'].astype(str) + "," + joinedlist_exon2['stop_exon2'].astype(str) + "," + joinedlist_exon3['stop_exon3'].astype(str) + "," + joinedlist_exon4['stop_exon4'].astype(str) + "," + joinedlist_exon5['stop_exon5'].astype(str) + "," + joinedlist_exon6['stop_exon6'].astype(str)+ "," + joinedlist_exon7['stop_exon7'].astype(str) + "," + joinedlist_exon8['stop_exon8'].astype(str) + "," + joinedlist_exon9['stop_exon9'].astype(str) + "," + joinedlist_exon10['stop_exon10'].astype(str) + ","
exon_startsdf = pd.DataFrame(exon_starts)
exon_starts_extract = exon_startsdf[0].str.split("None,", 1, expand = True)
exon_starts_extract.columns = ['exons_starts', 'throwaway']

exon_stopsdf = pd.DataFrame(exon_stops)
exon_stops_extract = exon_stopsdf[0].str.split("None,", 1, expand = True)
exon_stops_extract.columns = ['exons_stops', 'throwaway']

joinedlistfinal = joinedlist[0].astype(str) + "\t" + joinedlist[1].astype(str) + "\t" + joinedlist[2].astype(str) + "\t" + joinedlist[3].astype(str) + "\t" + joinedlist[4].astype(str) + "\t" + joinedlist[5].astype(str) + "\t" + joinedlist[6].astype(str) + "\t" + joinedlist[7].astype(str) + "\t" + joinedlist[8].astype(str) #+ " " + "exon_starts" + " \"" + exon_starts_extract['exons_starts'].astype(str) + "\"" + " " + "exon_stop" + " \"" + exon_stops_extract['exons_stops'].astype(str) + "\"" + ";"

fastafile = sys.argv[2]
#output_file = for fasta 
outfile1 = open(fastafile,'w')
for k in joinedlistfinal:
   outfile1.write(str(k) + "\n")
outfile1.close()

#formatting for genepred format
smorfid = joinedlist[8].str.split("\"", expand = True)


smorfid_genepred = smorfid[1]
smorfid_genepreddf = pd.DataFrame(smorfid_genepred)

exoncountdf = pd.DataFrame()
exoncountdf['count'] = exon_starts_extract['exons_starts'].str.count(",")

genepred_start = joinedlist[3].astype(int)-1
genepred_format = smorfid_genepreddf[1].astype(str) + "\t" + joinedlist[0].astype(str) + "\t" + joinedlist[6].astype(str) + "\t" + genepred_start.astype(str) + "\t" + joinedlist[4].astype(str) + "\t" + joinedlist[4].astype(str).astype(str) + "\t" + joinedlist[4].astype(str) + "\t" + exoncountdf['count'].astype(str) + "\t" + exon_starts_extract['exons_starts'].astype(str) + "\t" + exon_stops_extract['exons_stops'].astype(str)
genepred_format_list = genepred_format.values.tolist()

#genepred file output
genepred_file = sys.argv[3]
outfile_genepred = open(genepred_file,'w')
for j in genepred_format_list:
   outfile_genepred.write(str(j) + "\n")
outfile_genepred.close()
