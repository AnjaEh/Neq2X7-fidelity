#!/bin/bash

#Move through all data folders and apply cutadapt command to filter out the error enrichment site and trimm reads to flanking sequence of error enrichment site
for d in */; do
  echo "Processing directory: ${d}"
  if [ -f "${d}"*_R1_001.fastq.gz ]; then
    prefix="$(basename -- ${d})"
    prefix="${prefix:0:10}"
    cd "${d}"
    cutadapt -g file:../../adapter_read1.fasta -G file:../../adapter_read2.fasta --action=retain --no-indels -O=10 -o "${prefix}_R1_trimmed.fastq.gz" -p "${prefix}_R2_trimmed.fastq.gz" *_R1_001.fastq.gz *_R2_001.fastq.gz --untrimmed-output "${prefix}_R1_failed.fastq.gz" --untrimmed-paired-output "${prefix}_R2_failed.fastq.gz" --report=minimal
    cd ..
  fi
done


