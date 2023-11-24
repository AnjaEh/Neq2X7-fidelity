# Neq2X7-fidelity
Data analysis MagNIFI assay Neq2X7 DNA polymerase

Platform: Ubuntu 22.04.1 LTS  (Windows subsystem for Linux)  
Raw sequencing data found in folder: Neq2X7/FASTQ_Generation_2023-07-18 <br>
Well assignment: metadata.csv

Data analysis pipeline:  
1. Navigate to FASTQ_Generation_2023-07-18 directory
2. make scripts executable
3. excecute ./cutadapt.sh in activated cutadapt environment (cutadapt.yml)
4. exceute ./filter_align_results.sh in activated NGS-Linux2 environment (NGS-Linux2.yml)

* Results from analysis for Hernández-Rollán et al. (2023) are found in Alignment/pileup and Alignment/results.txt <br>

Hernández-Rollán et al. (2023) Neq2X7: a multi-purpose and open-source fusion DNA polymerase for advanced DNA engineering and diagnostics PCR
