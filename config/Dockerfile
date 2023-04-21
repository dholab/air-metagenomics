# Start with a base image
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install necessary packages
RUN apt-get update && \
    apt-get install -y python-is-python3 python3-pip r-base wget build-essential libz-dev libbz2-dev liblzma-dev \
    libncurses5-dev libncursesw5-dev libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev curl default-jre \
    samtools bcftools vcftools libhts-dev

# Install Miniconda
RUN curl -LO https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && \
    bash Miniconda3-latest-Linux-x86_64.sh -p /miniconda -b && \
    rm Miniconda3-latest-Linux-x86_64.sh

# Tinker with PATH to include SRA-tools
ENV PATH=/miniconda/bin:$PATH

# Create a new conda environment for sra-tools
RUN conda create -y --name sra_env -c defaults -c bioconda -c anaconda -c conda-forge -c r sra-tools=3.0.3 r-base r-tidyverse r-openxlsx

# Activate the sra_env environment and add sra-tools to the path
ENV PATH=/miniconda/envs/sra_env/bin:$PATH

# Install python packages
RUN pip3 install biopython

# Install cutadapt
RUN pip3 install cutadapt

# Make sure the python path is resolved despite using conda and pip
ENV PYTHONPATH="${PYTHONPATH}:/usr/local/lib/python3.10/site-packages"

# Install bbmap
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_39.01.tar.gz && \
    tar -xzf BBMap_39.01.tar.gz && \
    rm BBMap_39.01.tar.gz && \
    cd bbmap && \
    cp -r * /usr/local/bin

# Install minimap2
RUN wget https://github.com/lh3/minimap2/archive/v2.22.tar.gz && \
    tar -xzf v2.22.tar.gz && \
    rm v2.22.tar.gz && \
    cd minimap2-2.22 && \
    make && \
    cp minimap2 /usr/local/bin