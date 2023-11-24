# Neq2X7-fidelity
Data analysis MagNIFI assay Neq2X7 DNA polymerase

Platform: Ubuntu 22.04.1 LTS  (Windows subsystem for Linux)
Raw sequencing data found in folder: Neq2X7/FASTQ_Generation_2023-07-18  
Well assignment: metadata.csv

Data analysis pipeline:  
1. excecute FASTQ_Generation_2023-07-18/cutadapt.sh in activated cutadapt environment (cutadapt.yml)
2. exceute FASTQ_Generation_2023-07-18/filter_align_results.sh in activated NGS-Linux2 environment (NGS-Linux2.yml)

* Results from analysis for (Publication) are found in Alignment/pileup and Alignment/results.txt
