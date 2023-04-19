#!/usr/bin/env python3

import os
import sys
import gzip
from Bio import SeqIO
from multiprocessing import Pool, Lock

# Record the number of CPUs made available by nextflow
available_cpus = int(sys.argv[1])

# Set the directory containing the gzipped FASTQ files
fastq_dir = os.getcwd()

# Set the output directory for the converted FASTA files
fasta_dir = os.getcwd()

# Define a function to convert a single FASTQ file to FASTA
def convert_fastq_to_fasta(input_file):
    # Construct the path to the output file
    output_file = os.path.join(fasta_dir, os.path.basename(input_file).replace(".fastq.gz", ".fasta.gz"))

    # Open the input file
    with gzip.open(input_file, "rt") as in_handle:
        # Open the output file and lock it for writing
        with gzip.open(output_file, "wt") as out_handle, lock:
            # Use BioPython to convert the gzipped FASTQ file to FASTA format
            SeqIO.convert(in_handle, "fastq", out_handle, "fasta")

# Create a lock object to synchronize writing to the output file
lock = Lock()

# Loop through each gzipped FASTQ file in the directory
input_files = [os.path.realpath(filename) for filename in os.listdir(fastq_dir) if filename.endswith(".fastq.gz")]

# Process the files in parallel using a thread pool
with Pool(processes=available_cpus) as pool:
    pool.map(convert_fastq_to_fasta, input_files)
