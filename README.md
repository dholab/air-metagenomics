# Pathogen Agnostic Sequence Analysis

## Overview

NextFlow pipeline for quality controlling and aligning air cartridge sequence reads to a variety of possible pathogens. More narrative description coming soon.

## Quick Start

The two core software requirements to run this workflow are Nextflow to manage the data flow and Docker to provide software (though apptainer also works). If Nextflow and Docker are already installed on your system, no `git clone` commands or other setup is required. Simply run the workflow and reproduce our findings with:

```
nextflow run dholab/pathogen-agnostic-sequence-analysis -latest --samplesheet resources/samplesheet_27269.csv
```

With this command, Nextflow will automatically pull the workflow bundle of files from this GitHub repository and run it.

## Detailed Instructions

To run this workflow, start by downloading the workflow bundle to a working directory of your choice with `git clone`, like so:

```
git clone https://github.com/dholab/pathogen-agnostic-sequence-analysis.git .
```

When the workflow bundle has downloaded, you may need to set the workflow scripts to executable by running `chmod +x bin/*` in the command line.

If you haven't already, you will also need to install the Docker engine. The workflow pulls all the software it needs automatically from Docker Hub, which means you will never need to permanently install that software on your system. To install Docker, simply visit the Docker installation page at [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/). If you are running the workflow on a high-powered computing cluster where you do not have root permissions, you may also run the workflow with [Apptainer (formerly called Singularity)](https://apptainer.org/).

### Nextflow Installation

This workflow uses the [NextFlow](https://www.nextflow.io/) workflow manager. We recommend you install NextFlow to your system in one of the two following ways:

#### 1) Installation with Conda

1. Install the miniconda python distribution, if you haven't already: [https://docs.conda.io/en/latest/miniconda.html](https://docs.conda.io/en/latest/miniconda.html)
2. Install the `mamba` package installation tool in the command line:
   `conda install -y -c conda-forge mamba`
3. Install Nextflow to your base environment:
   `mamba install -c bioconda nextflow `

#### 2) Installation with curl

1. Run the following line in a directory where you'd like to install NextFlow, and run the following line of code:
   `curl -fsSL https://get.nextflow.io | bash`
2. Add this directory to your $PATH. If on MacOS, a helpful guide can be viewed [here](https://www.architectryan.com/2012/10/02/add-to-the-path-on-mac-os-x-mountain-lion/).

To double check that the installation was successful, type `nextflow -v` into the terminal. If it returns something like `nextflow version 22.10.0.5826`, you are set and ready to proceed. We developed the workflow with nextflow version 22.10.0.5826.

You are now ready to run the workflow!

### Running and Managing the Workflow

A core aspect of the workflow's logic is that it runs once for each sequencing run, with the requirement that each run contains at least one negative control. As such, to get the results presented in Ramuta et al. 2023, we ran the workflow twice with the following commands:

#### Sequencing run 1:

```
nextflow run main.nf --samplesheet resources/samplesheet_27269.csv
```

#### Sequencing run 2:

```
nextflow run main.nf --samplesheet resources/samplesheet_28210.csv
```
