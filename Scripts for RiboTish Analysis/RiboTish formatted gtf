#Before running this script, use the UCSC tool genePredtoGtf to create a gtf file from the genepred file ran in a previous script. Then this script will take the outputted gtf and format it for downstream filtering.
#Usage: python RiboTish_adjustcoordinates.py RiboTish_genepred_output.gtf RiboTish_formattedgtf.gtf

import sys
import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
inputgtf  = sys.argv[1]

genepred_gtf = pd.read_csv(inputgtf, sep="\t", header=None)
genepred_gtf.columns = ['chr', 'RiboTish', 'type', 'start', 'stop', 'dot1', 'strand', 'dot2', 'gene_id']
genepred_gtf_gene_id = genepred_gtf['gene_id'].str.split("\"", expand=True)
geneid = genepred_gtf_gene_id[1]
geneid = geneid.to_frame()
geneid.columns = ['gene_id_name']
genepred_gtf = pd.concat([genepred_gtf, geneid], axis = 1)

genepred_gtf_exons = genepred_gtf[genepred_gtf['type'].str.contains('exon')]
genepred_gtf_transcript = genepred_gtf[genepred_gtf['type'].str.contains('transcript')]

positivestrand_genepred_exons = genepred_gtf_exons[genepred_gtf_exons.strand == '+']
pos_exon_adjust = positivestrand_genepred_exons.loc[positivestrand_genepred_exons.groupby(["gene_id_name"])['stop'].idxmax()]
pos_exon_adjust['stop'] = pos_exon_adjust['stop'] - 3 
positivestrand_genepred_exons.update(pos_exon_adjust)
genepred_gtf_exons.update(positivestrand_genepred_exons)
#print(genepred_gtf_exons[:20])


negativestrand_genepred_exons = genepred_gtf_exons[genepred_gtf_exons.strand == '-']
neg_exon_adjust = negativestrand_genepred_exons.loc[negativestrand_genepred_exons.groupby(["gene_id_name"])['start'].idxmin()]
neg_exon_adjust['start'] = neg_exon_adjust['start'] + 3
negativestrand_genepred_exons.update(neg_exon_adjust)
genepred_gtf_exons.update(negativestrand_genepred_exons)
#print(neg_exon_adjust[:5])

frames = [genepred_gtf_exons, genepred_gtf_transcript]
adjusted_gtf = pd.concat(frames)
adjusted_gtf = pd.DataFrame(adjusted_gtf)


exceptions = adjusted_gtf.loc[(adjusted_gtf['start'] >= adjusted_gtf['stop'])]
exceptions = exceptions.assign(diff=exceptions['start'] - exceptions['stop'])

#if the start and stop are the same
exceptions['start'] = np.where(exceptions['diff'] == 0, exceptions['start'] - 1, exceptions['start'])

#if the start is 1 larger than stop
exceptions['start'] = np.where(exceptions['diff'] == 1, exceptions['start'] - 1, exceptions['start'])
exceptions['stop'] = np.where(exceptions['diff'] == 1, exceptions['stop'] + 1, exceptions['stop'])


#if the start is 2 larger than stop
exceptions['start'] = np.where(exceptions['diff'] == 2, exceptions['start'] - 3, exceptions['start'])
exceptions['stop'] = np.where(exceptions['diff'] == 2, exceptions['stop'] + 1, exceptions['stop'])

#if the start is 3 larger than stop
exceptions['start'] = np.where(exceptions['diff'] == 3, exceptions['start'] - 4, exceptions['start'])

#drop difference column in exceptions df and merge updated values to adjusted_gtf df
exceptions = exceptions.drop(columns=['diff'])
adjusted_gtf.update(exceptions)

adjusted_gtf = adjusted_gtf.drop(['gene_id_name'], axis=1)
adjusted_gtf['start'] = adjusted_gtf['start'].astype(int)
adjusted_gtf['stop'] = adjusted_gtf['stop'].astype(int)


adjust_gtf_final = adjusted_gtf['chr'].astype(str) + "\t" + adjusted_gtf['RiboTish'].astype(str).astype(str) + "\t" + adjusted_gtf['type'].astype(str) + "\t" + adjusted_gtf['start'].astype(str) + "\t" + adjusted_gtf['stop'].astype(str) + "\t" + adjusted_gtf['dot1'].astype(str) + "\t" + adjusted_gtf['strand'].astype(str) + "\t" + adjusted_gtf['dot2'].astype(str) + "\t" + adjusted_gtf['gene_id'].astype(str)
adjust_gtf_final = adjust_gtf_final.values.tolist()


#output_file formatted gtf
output_file = sys.argv[2]
outfile1 = open(output_file,'w')
for k in adjust_gtf_final:
   outfile1.write(str(k) + "\n")
outfile1.close()
