See our [wiki](https://github.com/gavinmdouglas/canola_pseudomonas_RNAseq/wiki) for a description of the code and workflows used for this project.

This repo contains the following folders and files:
* ```At_deseq2``` - DESeq2 output files after mapping read counts to Arabidopsis thaliana homologs.
* ```At_gene_sets``` - Gene sets expressed across specified days / tissues. Note that files that contain only a single day contain ALL genes expressed on that day, not just genes that didn't intersect with other days.
* ```blast_working``` - BLAST database files and output results for mapping Brassica napus genes to Arabidopsis thaliana homologs.
* ```mmquant_out``` - Output files after running mmquant to identify expression levels of Brassica napus genes (this program collapses duplicate genes into a single category)
* ```mmquant_out_AT``` - Output files after mapping mmquant results from B. napus to A. thaliana homologs.
* ```plots``` - Results and exploratory figures.
* ```scripts``` - R and Python scripts used for this project.
* ```tables``` - Result tables and mapping files produced for this project.
* ```canola_rnaseq_metadata_edit.txt``` - Metadata file for RNA-seq samples.
