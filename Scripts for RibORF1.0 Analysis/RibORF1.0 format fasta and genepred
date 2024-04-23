#input1 candidateORF.fa file produced from the ORFannotate function from RibORF1.0 is needed for input
#input2 from running RibORF1.0 is the repre.valid.pred.pvalue.parameters.txt file for collapsed ORF candidates was used for the analysis.
#input3 is the repre.valid.ORF.genepred.txt file  
#purpose of this script is to take the output of smORF candidates from RibORF1.0 and create a fasta and genepred file 
#Usage is python RibORF1.0_script.py candidateORF.fa repre.valid.pred.pvalue.parameters.txt repre.valid.ORF.genepred.txt fasta_formatted.fa repre.valid.ORF.genepred_formatted.txt

import sys
import numpy as np
import pandas as pd
inputfile  = sys.argv[1]
inputfile2 = sys.argv[2]
inputfile3 = sys.argv[3]
out = sys.argv[4]
out1 = sys.argv[5]

#candidate.fa
fasta = []
with open(inputfile,"r") as lines:
    for i in lines:
        if i.startswith('>'):
            newline = i.replace('\n','_')
            fasta.append(newline)
        else:    
            i1 = i.strip('')
            i2 = i1.replace('\n','')
            fasta.append(i2)


joinedfasta = ''.join(fasta)
splitjoinedfasta = joinedfasta.split(">")         
del splitjoinedfasta[0]

fasta_dataframe = pd.DataFrame(splitjoinedfasta, columns = ['smorf'])
fasta_dataframe = fasta_dataframe[~fasta_dataframe['smorf'].str.contains("PAR_")]
fasta_dataframe_split = fasta_dataframe['smorf'].str.split("_", expand = True)
fasta_dataframe_split.columns = ['orfID', 'AAsequence']




#input file from RibORF1.0
inputfile_inputpred = pd.read_csv(inputfile2, sep='\t', header='infer')
inputfile_inputpred.columns = ['orfID', 'chrom', 'strand', 'codon5', 'codon3', 'length', 'readNum', 'f1', 'f2', 'f3', 'entropy', 'Maxentropy', 'PME', 'codonNum', 'f1max', 'pred.pvalue']
#print(len(inputfile_inputpred))

inputfile_inputpred1 = inputfile_inputpred[inputfile_inputpred['pred.pvalue']>.7]

inputfile_inputpred2 = inputfile_inputpred1[inputfile_inputpred1['length']<454]

inputfile_inputpred2_1 = inputfile_inputpred2[inputfile_inputpred2['length']>18]

inputfile_inputpred3 = inputfile_inputpred2_1[inputfile_inputpred2_1['readNum']>10]

inputfile_inputpred3['codon5'] = np.where(inputfile_inputpred3['strand'] == '+', inputfile_inputpred3['codon5'] +1, inputfile_inputpred3['codon5'])
inputfile_inputpred3['codon3'] = np.where(inputfile_inputpred3['strand'] == '+', inputfile_inputpred3['codon3'] -3, inputfile_inputpred3['codon3'])

#negative strand
inputfile_inputpred3['codon5'] = np.where(inputfile_inputpred3['strand'] == '-', inputfile_inputpred3['codon5'] +4, inputfile_inputpred3['codon5'])

Merged_orfID = inputfile_inputpred3['orfID'].str.split(":", expand=True) 
Merged_formatted = Merged_orfID[0].astype(str) + inputfile_inputpred3['strand'].astype(str) + Merged_orfID[1].astype(str) + ":" + inputfile_inputpred3['codon5'].astype(str) + "-" + inputfile_inputpred3['codon3'].astype(str)

Merged_formatted=pd.DataFrame(Merged_formatted)

inputfile_inputpred3['ORFid_format'] = Merged_formatted
finalinput = inputfile_inputpred3.drop(columns=['chrom', 'strand', 'codon5', 'codon3', 'length', 'readNum', 'f1', 'f2', 'f3', 'entropy', 'Maxentropy', 'PME', 'codonNum', 'f1max', 'pred.pvalue'])

Mergedfasta_input = pd.merge(fasta_dataframe_split, finalinput, how='inner', on = 'orfID' )

finalfasta = '>' + Mergedfasta_input['ORFid_format'].astype(str) + '\n' + Mergedfasta_input['AAsequence'].astype(str)
#fastafile output
output_file = open(out, 'w')
for j in finalfasta :
    output_file.write(str(j) + "\n")
output_file.close

#genepred_inputfile from RibORF1.0
genepred_inputfile = pd.read_csv(inputfile3, sep="\t", header=None)
genepred_inputfile.columns = ['orfID', 'chr', 'strand', 'trans_start', 'trans_stop', 'cds_start', 'cds_stop', 'exon', 'exon_start', 'exon_stop']

genepred_merge = pd.merge(finalinput, genepred_inputfile, how='inner', on='orfID')
genepred_final = genepred_merge['ORFid_format'].astype(str) + "\t" + genepred_merge['chr'].astype(str) + "\t" + genepred_merge['strand'].astype(str) + "\t" + genepred_merge['trans_start'].astype(str) + "\t" + genepred_merge['trans_stop'].astype(str) + "\t" + genepred_merge['cds_start'].astype(str) + "\t" + genepred_merge['cds_stop'].astype(str) + "\t" + genepred_merge['exon'].astype(str) + "\t" + genepred_merge['exon_start'].astype(str) + "\t" + genepred_merge['exon_stop'].astype(str)
#out genepred file formatted
print(len(Mergedfasta_input))
print(len(genepred_merge))

output_file1 = open(out1, 'w')
for l in genepred_final :
    output_file1.write(str(l) + "\n")
output_file1.close
