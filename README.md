# Ribo-seq-Tool-Comparison-Scripts
Gregory Tong (4/23/24)

**Introduction**\
A list of scripts for formatting and comparing small open reading frame (smORF) prediction pipelines (RibORF0.1, RiboCode, RiboTish, ORFQuant, and RibORF1.0). To recreate the analysis, you must have the output files after running each of the ORF prediction tools and to run the following scripts in python3.

**List Comparison Scripts**\
In the list comparison scripts folder, there are two scripts that were used for formatting the fasta file to directly compare smORF results across the different analysis pipelines. Formatting_fasta.py is used for formatting the fasta files outputted from each of the tools excluding RibORF0.1. The other script is specifically for formatting the fastas generated from running RibORF0.1.

**RibORF0.1 Analysis**\
For producing smORF lists using RibORF0.1, details of the step by step command line code can be found [https://doi.org/10.1016/j.xpro.2023.102649] and [https://doi.org/10.1038/s41589-019-0425-0].

**RiboCode Analysis**\
Starting from the output files generated from running the RiboCode tool, the first step is to remove pseduoautosomal hits from the output files using sed 's/_PAR_Y//g'. Next, run the Formatting_RiboCode_gtf.py script using the RiboCode gtf output as the input to produce a formatted gtf. After, run the Generate_Fasta_RiboCode.py script using the RiboCode output.txt file as input to generate a formatted fasta. 

**Ribo-TISH Analysis**\
Starting from the output files generated from running the Ribo-TISH, the first script RiboTISH_genepred_conversion.py takes the RiboTish_output.txt file from RiboTish and produces a fasta file and genepred file of the ORF hits. The second script RiboTISH_formatted_fasta.py is used to convert the previously generated fasta file and convert to a formatted fasta file. Then using UCSC tools genepredtogtf, take the genepred file generated from the first script and output a gtf file. With that gtf file, run the script RiboTISH_formatted_gtf.py to produce a formatted gtf file.

**ORFquant**\
From running ORFquant, ORFquant produces an .Rdatafile that you can extract out with transcript ID name, seqnames, start, end, length, and strand and output into a csv file. Once that has been done, then you can run the first script ORFquant_formatted_fasta.py using the csv file as input and fasta that ORFquant outputs to generate a filtered fasta file. A blastp filter was used to filter for unannotated smORFs. Using the resulting blastp filtered fasta, run the second script ORFquant_formatted_gtf.py with the blastfiltered fasta, csv file containing (ORFid,  psites, ORF_pct_p_sites, ORF_pct_P_sites_pN, ORFs_pM, protein) from the .R data object and ORFquant produced gtf to generate a formatted gtf output file. 

**RibORF1.0**\
After running RibORF1.0, run the python script RibORF1.0_format_fasta_and_genepred.py using the candidateORF.fa file, collapsed repre.valid.pred.pvalue.parameters.txt, and the repre.valid.ORF.genepred.txt as input files. Otuput will be a formatted fasta file and genepred file. UCSC tools genePredToGtf function can take the genepred file and convert to a gtf file. Then the second script, RibORf1.0_format_gtf.py can be run with the resulting gtf file to create the resulting formatted gtf.

**Downstream Filtering**\
After running these scripts, additional filtering steps can be applied to the smORF lissts to remove ones that are annotated or overlap with known coding genes. Bedtools and NCBI blast can be used to filter smORFs that are overlapping coding regions or smORFs with similar amino acid sequences of known proteins. 


**Scripts For Analysis**\
Fraction_In_Frame_Reads.py script used to pull fraction 1 in frame reads from RibORF0.1 results to replicate analysis in the paper. The Long Read Analysis Commands folder contains a list of command line code used to process the long read data and generate transcriptome assemblies using StringTie. GFFcompare commands can also be found for generating transcriptome assembly statistics for long read and Illumina assemblies. TPM analysis contains the HOMER analyzeRepeat command to quantify coverage across smORF coding regions.



