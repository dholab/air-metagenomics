#!/usr/bin/env python3

import os
import sys
import gzip
from Bio import SeqIO

# Set the directory containing the gzipped FASTQ files
fastq_dir = sys.argv[1]

# Set the output directory for the converted FASTA files
fasta_dir = os.getcwd()

# Loop through each gzipped FASTQ file in the directory
for filename in os.listdir(fastq_dir):
    if filename.endswith(".fastq.gz"):
        # Construct the paths to the input and output files
        input_file = os.path.join(fastq_dir, filename)
        output_file = os.path.join(fasta_dir, filename.replace(".fastq.gz", ".fasta"))

        # Open the input and output files
        with gzip.open(input_file, "rt") as in_handle, open(
            output_file, "w"
        ) as out_handle:
            # Use BioPython to convert the gzipped FASTQ file to FASTA format
            SeqIO.convert(in_handle, "fastq", out_handle, "fasta")
