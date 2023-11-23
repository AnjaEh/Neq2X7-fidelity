#!/bin/bash

#Move through all data directories and apply paired read match check, general QC and alignment using bowtie2
for d in */; do
  echo "Processing directory: ${d}"
  if [ -f "${d}"*_R1_001.fastq.gz ]; then
    prefix="$(basename -- ${d})"
    prefix="${prefix:0:10}"
    cd "${d}"
    #Apply script to check for exact reverse complement match with paired read
    python ../paired_match_v2.py *_R1_trimmed.fastq.gz *_R2_trimmed.fastq.gz
    #Apply fastp program for final quality and length limit filtering
    fastp --in1 *_R1_trimmed_pass.fastq.gz --out1 "${prefix}_R1_filtered.fastq.gz" --in2 *_R2_trimmed_pass.fastq.gz --out2 "${prefix}_R2_filtered.fastq.gz" -A --unqualified_percent_limit=10 --length_required=50 --length_limit=70 --html "${prefix}_fastp.html"
    #Perform alignment with final, clean set of reads
    bowtie2 -x ../../templates_read1 -U *_R1_filtered.fastq.gz -S ../../Alignment/"${prefix}.sam" --local
    cd ..
  fi
done

#Move to alignment folder to process alignments
cd ../Alignment/

#Convert, sort, and index all .sam files in this subdirectory
for sam_file in *.sam; do
  # Get the name of the .sam file without the extension
  name="${sam_file%.*}"
  # Convert the .sam file to .bam format
  samtools view -bS "${sam_file}" | cat > "${name}.bam"
  # Sort the .bam file
  samtools sort "${name}.bam" > "${name}.sorted.bam"
  # Index the sorted .bam file
  samtools index "${name}.sorted.bam"
  #Mutation pileup and extract information about mismatches at position 45
  samtools mpileup "${name}.sorted.bam" -o ./pileup/"${name}_pileup.txt" -f templates_read1.fasta

  #Extract well identifier
  well="${name:7:10}"

  # Define the list of 4 identifiers
  identifiers=("T_a-rareT" "T_c-rareG" "T_t-rareA" "T_g-rareC")

  # Loop over each identifier and apply the commands
  # use samtools mpileup to extract data on matches and coverage at error enrichment site 
  for identifier in "${identifiers[@]}"; do
    samtools mpileup -f templates_read1.fasta -r "${identifier}:45-45" "${name}.sorted.bam" | cut -f 5 | grep -o '\.' | wc -l > match.txt
    samtools mpileup -f templates_read1.fasta -r "${identifier}:45-45" "${name}.sorted.bam" | cut -f 4 > depth.txt
    #samtools depth -r "${identifier}:30-60" "${name}.sorted.bam" | awk '{sum += $3} END {print sum/NR}' > depth.txt
    #write results line
    echo "$well" > well.txt
    echo "$identifier" > identifier.txt
    #Combine results line, append to full results table
    paste well.txt identifier.txt depth.txt match.txt > stats.txt
    cat stats.txt >> results.txt
  done
done
