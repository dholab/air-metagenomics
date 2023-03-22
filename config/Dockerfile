# Start with a base image
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND noninteractive
ENV TZ America/New_York

# Install necessary packages
RUN apt-get update && \
    apt-get install -y python3 python3-pip wget build-essential libz-dev libbz2-dev liblzma-dev \
    libncurses5-dev libncursesw5-dev libcurl4-openssl-dev libssl-dev zlib1g-dev default-jre \
    samtools bcftools vcftools libhts-dev

# Install sra-tools
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-ubuntu64.tar.gz && \
    tar -xzf sratoolkit.2.11.0-ubuntu64.tar.gz && \
    rm sratoolkit.2.11.0-ubuntu64.tar.gz && \
    cd sratoolkit.2.11.0-ubuntu64/bin && \
    cp fasterq-dump /usr/local/bin

# Install cutadapt
RUN pip3 install cutadapt

# Install bbmap
RUN wget https://downloads.sourceforge.net/project/bbmap/BBMap_38.90.tar.gz && \
    tar -xzf BBMap_38.90.tar.gz && \
    rm BBMap_38.90.tar.gz && \
    cd bbmap && \
    cp -r * /usr/local/bin

# Install minimap2
RUN wget https://github.com/lh3/minimap2/archive/v2.22.tar.gz && \
    tar -xzf v2.22.tar.gz && \
    rm v2.22.tar.gz && \
    cd minimap2-2.22 && \
    make && \
    cp minimap2 /usr/local/bin
