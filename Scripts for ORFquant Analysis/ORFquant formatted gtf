#Prerequisite: pull columns from ORFquant produced Rdata file contianing ORFid,  psites, ORF_pct_p_sites, ORF_pct_P_sites_pN, ORFs_pM, protein and create a csv file with those columns. 
#this script will take the previous outputted fasta, and gtf outputted by ORFquant and generate a formatted gtf with genomic relative coordinates to be used for downstream filtering.
#Usage: python ORFquant_gtf_format.py ORFquant_fasta_blastfiltered.fasta ORFquant_riboseq_info.csv ORFquant_for_gtf_output.gtf Formatted_ORFquantgtf.gtf

import sys
import pandas as pd
import numpy as np


smorfid = []
aaseq = []
#previous blastp filtered fasta list 
with open('pred.pvalue.parameters.unique.notInBlast.fasta') as ORFquant_Fasta_novelsmorfs:
    fasta_lines = ORFquant_Fasta_novelsmorfs.readlines()
    for k in fasta_lines:
        if k.startswith('>'):
            smorfid.append(k)
        else:
            aaseq.append(k)

Fastaformat_df = pd.DataFrame()
Fastaformat_df['smorfid'] = smorfid
Fastaformat_df['AAsequence'] = aaseq

Fastaformat_df['smorfid'] = Fastaformat_df['smorfid'].str.replace(">","")
Fastaformat_df['smorfid'] = Fastaformat_df['smorfid'].str.replace("\n","")
Fastaformat_df['AAsequence'] = Fastaformat_df['AAsequence'].str.replace("\n","")

Fastaformat_df = Fastaformat_df.loc[Fastaformat_df['AAsequence'].str.len() < 151]

chr_id = Fastaformat_df['smorfid']
chr_id = chr_id.str.replace("+", "-")
chr_id = chr_id.str.split("-", expand = True)
chr_extract = chr_id[1].str.split(":", expand = True)
fasta_chr = chr_extract[0]
fasta_chr_format = pd.DataFrame(fasta_chr)
fasta_chr_format.columns = ['chr']


fastamergedframes = [Fastaformat_df, fasta_chr_format]
fastaformat_mergeddf = pd.concat(fastamergedframes, axis = 1)



HeLaHiRes_orfinfo = pd.read_csv("ORFinfo.csv", sep = ",", header='infer')

Psite_adjust=HeLaHiRes_orfinfo['P_sites']*100
Psiteadjust_pd=pd.DataFrame(Psite_adjust)
Psiteadjust_pd=Psiteadjust_pd['P_sites'].astype(str)
Psiteadjust_pd=Psiteadjust_pd.str.split('.', expand=True)
Psiteadjust_format = Psiteadjust_pd[0].astype(int)
Psiteadjust_format = Psiteadjust_format/100
Psiteadjust_format_pd = pd.DataFrame(Psiteadjust_format)
Psiteadjust_format_pd.columns = ['P_sites']
Psiteadjust_format_pd_list = Psiteadjust_format_pd['P_sites'].astype(str)
Psitelist= Psiteadjust_format_pd_list.tolist()

Psitefinal_pd = pd.DataFrame()
Psitefinal_pd['P_sites'] = Psitelist

ORF_pct_P_sites = HeLaHiRes_orfinfo['ORF_pct_P_sites']*100
ORF_pct_P_sitesadjust_pd=pd.DataFrame(ORF_pct_P_sites)
ORF_pct_P_sitesadjust_pd=ORF_pct_P_sitesadjust_pd['ORF_pct_P_sites'].astype(str)
ORF_pct_P_sitesadjust_pd=ORF_pct_P_sitesadjust_pd.str.split('.', expand=True)
ORF_pct_P_sitesadjust_format = ORF_pct_P_sitesadjust_pd[0].astype(int)
ORF_pct_P_sitesadjust_format = ORF_pct_P_sitesadjust_format/100
ORF_pct_P_sitesadjust_format_pd = pd.DataFrame(ORF_pct_P_sitesadjust_format)
ORF_pct_P_sitesadjust_format_pd.columns = ['ORF_pct_P_sites']
ORF_pct_P_sitesadjust_format_pd_list = ORF_pct_P_sitesadjust_format_pd['ORF_pct_P_sites'].astype(str)
ORF_pct_P_siteslist= ORF_pct_P_sitesadjust_format_pd_list.tolist()
ORF_pct_P_sites_final_pd = pd.DataFrame()
ORF_pct_P_sites_final_pd['ORF_pct_P_sites'] = ORF_pct_P_siteslist


ORF_pct_P_sites_pN_adjust=HeLaHiRes_orfinfo['ORF_pct_P_sites_pN']*100
ORF_pct_P_sites_pNadjust_pd=pd.DataFrame(ORF_pct_P_sites_pN_adjust)
ORF_pct_P_sites_pNadjust_pd=ORF_pct_P_sites_pNadjust_pd['ORF_pct_P_sites_pN'].astype(str)
ORF_pct_P_sites_pNadjust_pd=ORF_pct_P_sites_pNadjust_pd.str.split('.', expand=True)
ORF_pct_P_sites_pNadjust_format = ORF_pct_P_sites_pNadjust_pd[0].astype(int)
ORF_pct_P_sites_pNadjust_format = ORF_pct_P_sites_pNadjust_format/100
ORF_pct_P_sites_pNadjust_format_pd = pd.DataFrame(ORF_pct_P_sites_pNadjust_format)
ORF_pct_P_sites_pNadjust_format_pd.columns = ['ORF_pct_P_sites_pN']
ORF_pct_P_sites_pNadjust_format_pd_list = ORF_pct_P_sites_pNadjust_format_pd['ORF_pct_P_sites_pN'].astype(str)
ORF_pct_P_sites_pNlist= ORF_pct_P_sites_pNadjust_format_pd_list.tolist()
ORF_pct_P_sites_pNfinal_pd = pd.DataFrame()
ORF_pct_P_sites_pNfinal_pd['ORF_pct_P_sites_pN'] = ORF_pct_P_sites_pNlist


ORFs_pM_adjust=HeLaHiRes_orfinfo['ORFs_pM']*100
ORFs_pMadjust_pd=pd.DataFrame(ORFs_pM_adjust)
ORFs_pMadjust_pd=ORFs_pMadjust_pd['ORFs_pM'].astype(str)
ORFs_pMadjust_pd=ORFs_pMadjust_pd.str.split('.', expand=True)
ORFs_pMadjust_format = ORFs_pMadjust_pd[0].astype(int)
ORFs_pMadjust_format = ORFs_pMadjust_format/100
ORFs_pMadjust_format_pd = pd.DataFrame(ORFs_pMadjust_format)
ORFs_pMadjust_format_pd.columns = ['ORFs_pM']
ORFs_pMadjust_format_pd_list = ORFs_pMadjust_format_pd['ORFs_pM'].astype(str)
ORFs_pMlist= ORFs_pMadjust_format_pd_list.tolist()

ORFs_pMfinal_pd = pd.DataFrame()
ORFs_pMfinal_pd['ORFs_pM'] = ORFs_pMlist
#print(ORFs_pMfinal_pd)

Protein = pd.DataFrame()
Protein['AAsequence'] = HeLaHiRes_orfinfo['Protein']

ORF_id = pd.DataFrame()
ORF_id['ORF_id'] = HeLaHiRes_orfinfo['ORF_id']

dframes = [Psitefinal_pd, ORF_pct_P_sites_final_pd, ORF_pct_P_sites_pNfinal_pd, ORFs_pMfinal_pd, Protein, ORF_id]
Formatted_orfinfo = pd.concat(dframes, axis=1)


#read in the gtf that gets outputted from orfquant
ORFquant_gtf = pd.read_csv("for_ORFquant_Detected_ORFs.gtf", sep = "\t", skiprows = 3, header=None)
#print(ORFquant_gtf[:5])
#extract all the CDS lines in the orfquant gtf
ORFquant_CDS = ORFquant_gtf[ORFquant_gtf[2].str.contains('CDS')]
#split accessory information field and extract out the ORFid
ORF_info = ORFquant_CDS[8].str.split("\"", expand=True)

ORF_info_accessory = ORF_info[12].str.split(";",expand=True)


#recreate a dataframe, replacing the cds field with exon
ORFquant_format = ORFquant_CDS[0].astype(str) + "\t" + ORFquant_CDS[1].astype(str) + "\t" + "exon" + "\t" + ORFquant_CDS[3].astype(str) + "\t" + ORFquant_CDS[4].astype(str) + "\t" + ORFquant_CDS[5].astype(str) + "\t" + ORFquant_CDS[6].astype(str) + "\t" + ORFquant_CDS[7].astype(str) + "\t" + ORF_info[11].astype(str) + "\t" + ORF_info_accessory[1].astype(str) + "\t" + ORF_info_accessory[2].astype(str) + "\t" + ORF_info_accessory[3].astype(str) + "\t" + ORF_info_accessory[4].astype(str)  
ORFquant_format = ORFquant_format.values.tolist()
ORFquant_gtf_df = pd.DataFrame(x.split('\t') for x in ORFquant_format)
ORFquant_gtf_df.columns = ['chr', 'ORFquant', 'exon', 'start', 'stop', 'dot1', 'strand', 'dot2', 'ORF_id', 'P_sites', 'ORF_pct_P_sites', 'ORF_pct_P_sites_pN', 'ORFs_pM']
ORFquant_gtf_df['P_sites'] = ORFquant_gtf_df['P_sites'].str.replace("P_sites ", "")
ORFquant_gtf_df['ORF_pct_P_sites'] = ORFquant_gtf_df['ORF_pct_P_sites'].str.replace("ORF_pct_P_sites ", "")
ORFquant_gtf_df['ORF_pct_P_sites_pN'] = ORFquant_gtf_df['ORF_pct_P_sites_pN'].str.replace("ORF_pct_P_sites_pN ", "")
ORFquant_gtf_df['ORFs_pM'] = ORFquant_gtf_df['ORFs_pM'].str.replace("ORFs_pM ", "")

ORFquant_gtf_Psites_adjust = pd.Series([])
ORFquant_gtf_Psites_adjust=ORFquant_gtf_df['P_sites'].str.strip()
ORFquant_gtf_Psites_adjust=ORFquant_gtf_Psites_adjust.astype(float)*100
ORFquant_gtf_Psitesadjust_pd=pd.DataFrame(ORFquant_gtf_Psites_adjust)
ORFquant_gtf_Psitesadjust_pd=ORFquant_gtf_Psitesadjust_pd['P_sites'].astype(str)
ORFquant_gtf_Psitesadjust_pd=ORFquant_gtf_Psitesadjust_pd.str.split('.', expand=True)
ORFquant_gtf_Psitesadjust_format = ORFquant_gtf_Psitesadjust_pd[0].astype(int)
ORFquant_gtf_Psitesadjust_format = ORFquant_gtf_Psitesadjust_format/100
ORFquant_gtf_Psitesadjust_format_pd = pd.DataFrame(ORFquant_gtf_Psitesadjust_format)
ORFquant_gtf_Psitesadjust_format_pd.columns = ['P_sites']
ORFquant_gtf_Psitesadjust_format_pd_list = ORFquant_gtf_Psitesadjust_format_pd['P_sites'].astype(str)
ORFquant_gtf_Psiteslist= ORFquant_gtf_Psitesadjust_format_pd_list.tolist()

ORFquant_gtf_Psitesfinal_pd = pd.DataFrame()
ORFquant_gtf_Psitesfinal_pd['P_sites'] = ORFquant_gtf_Psiteslist


ORFquant_gtf_ORF_pct_P_sites_adjust = pd.Series([])
ORFquant_gtf_ORF_pct_P_sites_adjust=ORFquant_gtf_df['ORF_pct_P_sites'].str.strip()
ORFquant_gtf_ORF_pct_P_sites_adjust=ORFquant_gtf_ORF_pct_P_sites_adjust.astype(float)*100
ORFquant_gtf_ORF_pct_P_sitesadjust_pd=pd.DataFrame(ORFquant_gtf_ORF_pct_P_sites_adjust)
ORFquant_gtf_ORF_pct_P_sitesadjust_pd=ORFquant_gtf_ORF_pct_P_sitesadjust_pd['ORF_pct_P_sites'].astype(str)
ORFquant_gtf_ORF_pct_P_sitesadjust_pd=ORFquant_gtf_ORF_pct_P_sitesadjust_pd.str.split('.', expand=True)
ORFquant_gtf_ORF_pct_P_sitesadjust_format = ORFquant_gtf_ORF_pct_P_sitesadjust_pd[0].astype(int)
ORFquant_gtf_ORF_pct_P_sitesadjust_format = ORFquant_gtf_ORF_pct_P_sitesadjust_format/100
ORFquant_gtf_ORF_pct_P_sitesadjust_format_pd = pd.DataFrame(ORFquant_gtf_ORF_pct_P_sitesadjust_format)
ORFquant_gtf_ORF_pct_P_sitesadjust_format_pd.columns = ['ORF_pct_P_sites']
ORFquant_gtf_ORF_pct_P_sitesadjust_format_pd_list = ORFquant_gtf_ORF_pct_P_sitesadjust_format_pd['ORF_pct_P_sites'].astype(str)
ORFquant_gtf_ORF_pct_P_siteslist= ORFquant_gtf_ORF_pct_P_sitesadjust_format_pd_list.tolist()
ORFquant_gtf_ORF_pct_P_sitesfinal_pd = pd.DataFrame()
ORFquant_gtf_ORF_pct_P_sitesfinal_pd['ORF_pct_P_sites'] = ORFquant_gtf_ORF_pct_P_siteslist


ORFquant_gtf_ORF_pct_P_sites_pN_adjust = pd.Series([])
ORFquant_gtf_ORF_pct_P_sites_pN_adjust=ORFquant_gtf_df['ORF_pct_P_sites_pN'].str.strip()
ORFquant_gtf_ORF_pct_P_sites_pN_adjust=ORFquant_gtf_ORF_pct_P_sites_pN_adjust.astype(float)*100
ORFquant_gtf_ORF_pct_P_sites_pNadjust_pd=pd.DataFrame(ORFquant_gtf_ORF_pct_P_sites_pN_adjust)
ORFquant_gtf_ORF_pct_P_sites_pNadjust_pd=ORFquant_gtf_ORF_pct_P_sites_pNadjust_pd['ORF_pct_P_sites_pN'].astype(str)
ORFquant_gtf_ORF_pct_P_sites_pNadjust_pd=ORFquant_gtf_ORF_pct_P_sites_pNadjust_pd.str.split('.', expand=True)
ORFquant_gtf_ORF_pct_P_sites_pNadjust_format = ORFquant_gtf_ORF_pct_P_sites_pNadjust_pd[0].astype(int)
ORFquant_gtf_ORF_pct_P_sites_pNadjust_format = ORFquant_gtf_ORF_pct_P_sites_pNadjust_format/100
ORFquant_gtf_ORF_pct_P_sites_pNadjust_format_pd = pd.DataFrame(ORFquant_gtf_ORF_pct_P_sites_pNadjust_format)
ORFquant_gtf_ORF_pct_P_sites_pNadjust_format_pd.columns = ['ORF_pct_P_sites_pN']
ORFquant_gtf_ORF_pct_P_sites_pNadjust_format_pd_list = ORFquant_gtf_ORF_pct_P_sites_pNadjust_format_pd['ORF_pct_P_sites_pN'].astype(str)
ORFquant_gtf_ORF_pct_P_sites_pNlist= ORFquant_gtf_ORF_pct_P_sites_pNadjust_format_pd_list.tolist()
ORFquant_gtf_ORF_pct_P_sites_pNfinal_pd = pd.DataFrame()
ORFquant_gtf_ORF_pct_P_sites_pNfinal_pd['ORF_pct_P_sites_pN'] = ORFquant_gtf_ORF_pct_P_sites_pNlist


ORFquant_gtf_ORFs_pM_adjust = pd.Series([])
ORFquant_gtf_ORFs_pM_adjust=ORFquant_gtf_df['ORFs_pM'].str.strip()
ORFquant_gtf_ORFs_pM_adjust=ORFquant_gtf_ORFs_pM_adjust.astype(float)*100
ORFquant_gtf_ORFs_pMadjust_pd=pd.DataFrame(ORFquant_gtf_ORFs_pM_adjust)
ORFquant_gtf_ORFs_pMadjust_pd=ORFquant_gtf_ORFs_pMadjust_pd['ORFs_pM'].astype(str)
ORFquant_gtf_ORFs_pMadjust_pd=ORFquant_gtf_ORFs_pMadjust_pd.str.split('.', expand=True)
ORFquant_gtf_ORFs_pMadjust_format = ORFquant_gtf_ORFs_pMadjust_pd[0].astype(int)
ORFquant_gtf_ORFs_pMadjust_format = ORFquant_gtf_ORFs_pMadjust_format/100
ORFquant_gtf_ORFs_pMadjust_format_pd = pd.DataFrame(ORFquant_gtf_ORFs_pMadjust_format)
ORFquant_gtf_ORFs_pMadjust_format_pd.columns = ['ORFs_pM']
ORFquant_gtf_ORFs_pMadjust_format_pd_list = ORFquant_gtf_ORFs_pMadjust_format_pd['ORFs_pM'].astype(str)
ORFquant_gtf_ORFs_pMlist= ORFquant_gtf_ORFs_pMadjust_format_pd_list.tolist()

ORFquant_gtf_ORFs_pMfinal_pd = pd.DataFrame()
ORFquant_gtf_ORFs_pMfinal_pd['ORFs_pM'] = ORFquant_gtf_ORFs_pMlist


ORFquant_gtf_ORF_id = ORFquant_gtf_df['ORF_id']
ORFquant_gtf_ORF_id_pd = pd.DataFrame(ORFquant_gtf_ORF_id).astype(str)
ORFquant_gtf_ORF_id_pd.columns = ['ORF_id']



#create a gtf dataframe now with the appropriate formatted columns
ORFquant_gtf_df = ORFquant_gtf_df.drop(['P_sites', 'ORF_pct_P_sites', 'ORF_pct_P_sites_pN', 'ORFs_pM'], axis=1)
formatgtf_df = [ORFquant_gtf_df, ORFquant_gtf_Psitesfinal_pd, ORFquant_gtf_ORF_pct_P_sitesfinal_pd, ORFquant_gtf_ORF_pct_P_sites_pNfinal_pd, ORFquant_gtf_ORFs_pMfinal_pd]
formatted_orf_gtf = pd.concat(formatgtf_df, axis=1)

FastaMerged_df = pd.merge(Formatted_orfinfo, fastaformat_mergeddf, how='inner', on=['AAsequence'])
FastaMerged_df = FastaMerged_df.dropna(axis=0)


Final_df_gtf_withsmorfid = pd.merge(FastaMerged_df, formatted_orf_gtf, how='inner', on=['P_sites', 'ORF_pct_P_sites', 'ORF_pct_P_sites_pN', 'ORFs_pM', 'chr'])
Final_df_gtf_withsmorfid = Final_df_gtf_withsmorfid.dropna(axis=0)

ORFquant_outputgtf = Final_df_gtf_withsmorfid['chr'].astype(str) + "\t" + Final_df_gtf_withsmorfid['ORFquant'].astype(str) + "\t" + Final_df_gtf_withsmorfid['exon'].astype(str) + "\t" + Final_df_gtf_withsmorfid['start'].astype(str) + "\t" + Final_df_gtf_withsmorfid['stop'].astype(str) + "\t" + Final_df_gtf_withsmorfid['dot1'].astype(str) + "\t" + Final_df_gtf_withsmorfid['strand'].astype(str) + "\t" + Final_df_gtf_withsmorfid['dot2'].astype(str) + "\t" + "gene_id " + "\"" + Final_df_gtf_withsmorfid['smorfid'].astype(str) + "\"; " + "transcript_id " + "\"" + Final_df_gtf_withsmorfid['smorfid'].astype(str) + "\"" + ";"
ORFquant_outputgtf = ORFquant_outputgtf.tolist()

out = "ORFquant_gtfformatted_gtf.gtf"
output_file = open(out, 'w')
for i in ORFquant_outputgtf:
    output_file.write(str(i) + "\n")
output_file.close

