#!/usr/bin/env python3

import os
import sys
import gzip
from Bio import SeqIO

# Set the directory containing the gzipped FASTQ files
fastq_dir = os.getcwd()

# Set the output directory for the converted FASTA files
fasta_dir = os.getcwd()

# Loop through each gzipped FASTQ file in the directory
for filename in os.listdir(fastq_dir):
    if filename.endswith(".fastq.gz"):
        # Construct the paths to the input and output files
        input_file = os.path.realpath(filename)
        output_file = os.path.join(fasta_dir, filename.replace(".fastq.gz", ".fasta.gz"))

        # Open the input and output files
        with gzip.open(input_file, "rt") as in_handle, gzip.open(
            output_file, "wt"
        ) as out_handle:
            # Use BioPython to convert the gzipped FASTQ file to FASTA format
            SeqIO.convert(in_handle, "fastq", out_handle, "fasta")

# Close the input and output files
in_handle.close()
out_handle.close()