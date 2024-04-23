# Ribo-seq-Tool-Comparison-Scripts
Gregory Tong (4/23/24)

**Introduction**\
A list of scripts for formatting and comparing small open reading frame (smORF) prediction pipelines (RibORF0.1, RiboCode, RiboTish, ORFQuant, and RibORF1.0). To recreate the analysis, you must have the output files after running each of the ORF prediction tools and to run the following scripts in python3.

**List Comparison Scripts**
In the list comparison scripts folder, there are two scripts that were used for formatting the fasta file to directly compare smORF results across the different analysis pipelines. Formatting_fasta.py is used for formatting the fasta files outputted from each of the tools excluding RibORF0.1. The other script is specifically for formatting the fastas generated from running RibORF0.1.

**RibORF0.1 Analysis**\
For producing smORF lists using RibORF0.1, details of the step by step command line code can be found [https://doi.org/10.1016/j.xpro.2023.102649] and [https://doi.org/10.1038/s41589-019-0425-0].

**RiboCode Analysis**
Starting from the output files generated from running the RiboCode tool, the first step is to remove pseduoautosomal hits from the output files using sed 's/_PAR_Y//g'. Next, run the Formatting_RiboCode_gtf.py script using the RiboCode gtf output as the input to produce a formatted gtf. After, run the Generate_Fasta_RiboCode.py script using the RiboCode output.txt file as input to generate a formatted fasta. 

**RiboTish Analysis**








