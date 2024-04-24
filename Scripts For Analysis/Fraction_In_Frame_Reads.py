#Script for taking frame 1 fraction of in frame reads from RibORF0.1 smORF output list
#Input 1: Output files from RibORF0.1 (pred.pvalue.parameters.txt)
#Input 2: smORF list from RibORF0.1
#Out2: csv file with fraction of reads in frame 1 and smORF IDs
import sys
import pandas as pd
import numpy as np

input1 = sys.argv[1]
intpu2 = sys.argv[2]
out2 = sys.argv[3]

RibORF_output = pd.read_csv(input1, sep="\t", header='infer')
geneID = RibORF_output['geneID'].str.split("_", expand = True)
geneID_split1 = geneID[0].str.split("chr", expand=True)
strand = geneID_split1[0].str.strip().str[-1]

geneID_formatted = strand + geneID_split1[1]

RibORF_inputparameters = pd.DataFrame()
RibORF_inputparameters['geneID'] = geneID_formatted
RibORF_inputparameters['frac1'] = RibORF_output['f1']

smORFList = pd.read_csv(input2, header=None)
smORFList_split = smORFList[0].str.split("_", expand=True)

smORFList_df = pd.DataFrame()
smORFList_df['geneID'] = smORFList_split[0]

merged_dataframe = pd.merge(RibORF_inputparameters, smORFList_df, how='inner', on=['geneID'])
uniqueID = merged_dataframe.drop_duplicates(subset = ['geneID'])

uniqueID.to_csv(out2, index=False)
