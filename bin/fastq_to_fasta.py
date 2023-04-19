#!/usr/bin/env python3

import os
import sys
import gzip
from Bio import SeqIO
from multiprocessing import Pool, Lock, Manager
from functools import partial

# Record the number of CPUs made available by nextflow
available_cpus = int(sys.argv[1])

# Define the number of reads to process before writing to file
batch_size = int(sys.argv[2])

# Set the directory containing the gzipped FASTQ files
fastq_dir = os.getcwd()

# Set the output directory for the converted FASTA files
fasta_dir = os.getcwd()

# Define a function that accepts a file writing lock and a batch of input files
def convert_fastq_to_fasta(lock, input_files):
    for input_file in input_files:
        
        # Construct the path to the output file
        output_file = os.path.join(fasta_dir, os.path.basename(input_file).replace(".fastq.gz", ".fasta.gz"))

        # Open the input and output files and the lock
        with gzip.open(input_file, "rt") as in_handle, gzip.open(output_file, "wt") as out_handle, lock:
            # Use BioPython to convert the gzipped FASTQ file to FASTA format
            SeqIO.convert(in_handle, "fastq", out_handle, "fasta")

# Define a function that handles multithreading and file locking for each batch of input files
def process_input_files(input_files, num_threads, batch_size):
    
    # Create a file lock for this thread
    lock = manager.Lock()

    # Split the input_files into chunks of 500
    input_files_chunks = [input_files[i:i+batch_size] for i in range(0, len(input_files), batch_size)]

    # Create a pool of worker processes
    with Pool(processes=num_threads) as pool:
        # Use the pool to map the conversion function to each chunk of input files
        pool.map(partial(convert_fastq_to_fasta, lock), input_files_chunks)

if __name__ == '__main__':
    # Get a list of all the gzipped FASTQ files in the directory
    input_files = [os.path.realpath(filename) for filename in os.listdir(fastq_dir) if filename.endswith(".fastq.gz")]
    
    # Create a manager object to manage the locks across different threads
    with Manager() as manager:
        # Process the input files in parallel using a thread pool
        process_input_files(input_files, num_threads=available_cpus, batch_size=batch_size)
