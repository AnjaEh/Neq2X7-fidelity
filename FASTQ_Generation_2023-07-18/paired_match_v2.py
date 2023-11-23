import sys
import gzip
from Bio import SeqIO

def reverse_complement(seq):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}
    return ''.join([complement[base] for base in reversed(seq)])

def check_reverse_complement(read1, read2):
    """
    Check if read1 is the reverse complement of read2, and vice versa.
    Returns True if they are, False otherwise.
    """
    # Check if the reads are the same length
    if len(read1) != len(read2):
        global length_mismatch_count
        length_mismatch_count += 1
        return False
    else:
        # Check if read1 is the reverse complement of read2
        if read1.reverse_complement() == read2:
            return True
        # Check if read2 is the reverse complement of read1
        elif read2.reverse_complement() == read1:
            return True
        else:
            return False

if len(sys.argv) != 3:
    print("Usage: python script.py read1.fastq.gz read2.fastq.gz")
    sys.exit(1)

# Use the command-line arguments as input file names
read1_file = sys.argv[1]
read2_file = sys.argv[2]

output_file_r1_pass = f"{read1_file[:-9]}_pass.fastq.gz"
output_file_r2_pass = f"{read2_file[:-9]}_pass.fastq.gz"
output_file_r1_fail = f"{read1_file[:-9]}_fail.fastq.gz"
output_file_r2_fail = f"{read2_file[:-9]}_fail.fastq.gz"

length_mismatch_count = 0
pass_count = 0
fail_count = 0


# Open the input files with gzip.open to handle compressed input files
with gzip.open(read1_file, "rt") as handle1, gzip.open(read2_file, "rt") as handle2, gzip.open(output_file_r1_pass, "wt") as out_handle_r1_pass, gzip.open(output_file_r2_pass, "wt") as out_handle_r2_pass, gzip.open(output_file_r1_fail, "wt") as out_handle_r1_fail, gzip.open(output_file_r2_fail, "wt") as out_handle_r2_fail:
    #num_exact_match_reads = 0
    for rec1, rec2 in zip(SeqIO.parse(handle1, "fastq"), SeqIO.parse(handle2, "fastq")):
        # Check if the reads are reverse complements of each other
        if check_reverse_complement(rec1.seq, rec2.seq):
            # Write the reads to the output file for passing reads
            SeqIO.write(rec1, out_handle_r1_pass, "fastq")
            SeqIO.write(rec2, out_handle_r2_pass, "fastq")
            pass_count += 1
        else:
            # Write the reads to the output file for failing reads
            SeqIO.write(rec1, out_handle_r1_fail, "fastq")
            SeqIO.write(rec2, out_handle_r2_fail, "fastq")
            fail_count += 1

print(f"{pass_count} read pairs passed the reverse complement filter and were written to {output_file_r1_pass} and the paired file.")
print(f"{fail_count} read pairs failed the reverse complement filter and were written to {output_file_r1_fail} and the paired file.")
print(f"{length_mismatch_count} read pairs were not the same length.")