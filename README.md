# Pathogen Agnostic Sequence Analysis

Developed for [Minor and Ramuta et al. 2023, _Metagenomic sequencing detects human respiratory and enteric viruses in air samples collected from congregate settings_](https://www.nature.com/articles/s41598-023-48352-6)

## Overview

Pathogen surveillence data from the COVID-19 pandemic overwhelmingly came from individual-testing, SARS-CoV-2 specific assays, such as nasal swab qPCR or rapid antigen tests. However, as individual testing becomes less common in 2023, researchers and public health officials are pivoting toward more population-level, environmental sample-based surveillance. Even still, the laboratory preparations and bioinformatics developed for environmental samples in the COVID-19 pandemic have tended to be pathogen-specific, leaving us blind to the majority of circulating pathogens that may or may not be present in an environmental sample.

In _Ramuta et al. 2023_, we describe a pathogen-agnostic, metagenomic approach to detecting any pathogen in the air. Using genetic material from air samplers distributed across the Upper Midwest, we show that air samplers can detect a variety of respiratory and enteric pathogens through time, without enriching for any one pathogen in particular. This approach could inform public health decision-making with data on all circulating pathogens, not just those that are expected.

The goal of the Nextflow pipeline and associated files in this repository is to make our results reproducible, and to make the methods we used to generate those results as portable and transparent as possible. We invite anyone interested to view our open data portal at [go.wisc.edu/h8mn47](https://go.wisc.edu/h8mn47), and to make use of our data in [NCBI BioProject # PRJNA950127](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA950127).

## Quick Start

The two core software requirements to run this workflow are Nextflow to manage the data flow and Docker to provide software (though Apptainer also works). If Nextflow and Docker are already installed on your system, no `git clone` commands or other setup is required. Simply run the workflow and reproduce our findings with:

```
nextflow run dholab/pathogen-agnostic-sequence-analysis -latest --samplesheet resources/samplesheet_27269.csv
```

With this command, Nextflow will automatically pull the workflow bundle of files from this GitHub repository and run it.

## Table of Contents

- [Detailed Instructions](https://github.com/dholab/pathogen-agnostic-sequence-analysis#detailed-instructions)
  - [Nextflow Installation](https://github.com/dholab/pathogen-agnostic-sequence-analysis#nextflow-installation)
    - [Installation with Conda](https://github.com/dholab/pathogen-agnostic-sequence-analysis#1-installation-with-conda)
    - [Installation with curl](https://github.com/dholab/pathogen-agnostic-sequence-analysis#2-installation-with-curl)
  - [Running and Managing the Workflow](https://github.com/dholab/pathogen-agnostic-sequence-analysis#running-and-managing-the-workflow)
- [Workflow Configuration](https://github.com/dholab/pathogen-agnostic-sequence-analysis#workflow-configuration)
- [Workflow Steps](https://github.com/dholab/pathogen-agnostic-sequence-analysis#workflow-configuration)
- [Acknowledgements](https://github.com/dholab/pathogen-agnostic-sequence-analysis#acknowledgements)
- [Citing the Workflow](https://github.com/dholab/pathogen-agnostic-sequence-analysis#citing-the-workflow)

## Detailed Instructions

To run this workflow, start by downloading the workflow bundle to a working directory of your choice with `git clone`, like so:

```
git clone https://github.com/dholab/pathogen-agnostic-sequence-analysis.git .
```

When the workflow bundle has downloaded, you may need to set the workflow scripts to executable by running `chmod +x bin/*` in the command line.

If you haven't already, you will also need to install the Docker engine. The workflow pulls all the software it needs automatically from Docker Hub, which means you will never need to permanently install that software on your system. To install Docker, simply visit the Docker installation page at [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/). If you are running the workflow on a high-powered computing (HPC) cluster where you do not have root permissions, you may also run the workflow with [Apptainer (formerly called Singularity)](https://apptainer.org/).

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

To run properly, the workflow requires at least 16 gigabytes of RAM, with as much as 250 gigabytes of available disk space. Your mileage may vary with less available memory and disk space.

A core aspect of the workflow's logic is that it runs once for each sequencing run, with the requirement that each run contains at least one negative control. As such, to get the results presented in Ramuta et al. 2023, we ran the workflow twice with the following commands:

#### Sequencing run 1:

```
nextflow run main.nf --samplesheet resources/samplesheet_27269.csv
```

#### Sequencing run 2:

```
nextflow run main.nf --samplesheet resources/samplesheet_28210.csv
```

As with any Nextflow workflow, you may resume an interrupted workflow run with the `-resume` command line argument, like so:

```
nextflow run main.nf --samplesheet resources/samplesheet_28210.csv -resume
```

You may also run the workflow in the background with the `-bg` flag:

```
nextflow -bg run main.nf --samplesheet resources/samplesheet_28210.csv
```

## Workflow Configuration

The most important file to configure this workflow for your own purposes is the sample sheet. This CSV-formatted table has three columns:

- `raw_read_label`: This is either a) a sample identifier that is part of the file names for locally stored FASTQs, or b) an SRA accession. If it is the former, the workflow will merge all FASTQ files that have this label. If it is an SRA accession, it will automatically pull the FASTQ from SRA. We recommend using the latter functionality and supporting open sequence databases.
- `sample_id`: An identifier you'd like to use for each sample. This is generally a more readable label, but could be whatever you want.
- `parent_dir`: The parent directory file path for each sample's files. This column corresponds to using a `raw_read_label` for a locally stored file. If the `raw_read_label` column is all SRA accessions, you may leave this column blank.

Additionally, there are a number of useful settings built into this workflow. The first is a `low_disk_mode`, which prevents the workflow from copying large output files into the results directory. This will roughly half the amount of disk space used by the workflow. It is invoked as a simple true/false boolean, like so:

```
nextflow run main.nf --samplesheet resources/samplesheet_28210.csv --low_disk_mode true
```

The workflow also allows you to specify an adapter sequence with the parameter `adapter_seq`, a reference sequence file with the parameter `pathogen_ref`, and a results directory with `results`, as demonstrated in the more complex command below:

```
nextflow run main.nf \
--samplesheet resources/samplesheet.csv \
--pathogen_ref /absolute/path/to/reference.fasta \
--results ~/Downloads/new_results \
--adapter_seq TTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT \
--low_disk_mode true
```

Finally, we've made the workflow sensitive to running in a variety of different compute environments. One way we do this is with the `--max_local_cpus` parameter, which you can use to constrain the total number of CPUs that the workflow draws from. Note, however, that using fewer CPUs will mean a longer runtime.

If you are running on an HPC cluster, we recommend you use the `hpc_cluster` profile, like so:

```
nextflow run main.nf \
-profile hpc_cluster \
--samplesheet relative/path/on/node/to/samplesheet.csv \
--pathogen_ref relative/path/on/node/to/reference.fasta \
--results new_results \
--adapter_seq TTTTTTTTTTAATGTACTTCGTTCAGTTACGTATTGCT \
--low_disk_mode true
```

## Workflow Steps

Broadly speaking, the workflow goes through three phases: 1) Sequence read collection and QC, 2) Sequence read filtering, and 3) Read alignment and reporting. All steps are visualized in [the workflow visualization](https://github.com/dholab/pathogen-agnostic-sequence-analysis/blob/main/workflow-visualization.png) in this repository.

In phase one, the workflow goes through the following steps:

- `FIND_AND_MERGE_FASTQS`: Here, the workflow either finds demultiplexed reads locally and merges them into one file per sample, or pulls reads automatically from SRA.
- `SAMPLE_QC`: Next, reads are filtered to a minimum length and quality score, and adapters are trimmed from either end.
- `FIND_NTC`: Negative control reads, which must be labeled with "NTC\_" in the sample ID, are identified and converted to the more lightweight FASTA format.
- `CONVERT_TO_FASTA`: Reads corresponding to samples, not controls, are converted to the more lightweight FASTA format.

In phase two, the data moves through these steps:

- `DOWNLOAD_CONTAMINANTS`: Five FASTA files, each of which contain exemplar sequences for a different category of common contaminants, are downloaded automatically from the Ramuta et al. 2023 open data sharing portal.
- `DECOMPRESS_CONTAMINANTS`: Contaminant FASTAs are downloaded in a compressed tarball to speed up the transfer. Here, they are decompressed to make them accessible downstream.
- `REMOVE_CONTAMINANTS`: Reads are then mapped to each contaminant FASTA. Only reads that do not map to contaminants—and are thus likely derived from genetic material sampled in the air cartridge—are retained downstream.
- `REMOVE_NTC`: Finally, reads are mapped to their sequencing run's negative control reads, which removed any contamination present from library preparation.

And in phase three, these steps finish off the workflow:

- `MAP_TO_REFSEQS`: Clean reads are then mapped to a FASTA file containing an arbitrary number of possible pathogen sequences.
- `RECORD_HITS`: The number of reads that successfully map to any reference sequence are recorded in a text file here.
- `MAKE_VIRUS_LOOKUP`: Here, the FASTA deflines from the pathogen reference file are parsed into NCBI RefSeq Accession and common name.
- `GENERATE_PIVOT_TABLE`: Finally, the workflow finishes by generating an Excel Pivot Table to report the positive hits, along with read support, in each sample.

## Acknowledgements

We thank the developers of Nextflow, BBTools, minimap2, and samtools for the critical subcomponents of this workflow. We also thank the students and staff at the AIDS Vaccine Research Laboratory of University of Wisconsin—Madison for their tireless work and helpful advice throughout the COVID-19 pandemic. And finally, we thank the entire ORegon CHild Absenteeism Due To Respiratory Disease Study (ORCHARDS) team, particularly Jonathan Temte, for making this study possible, and more broadly, for supporting air sampling as an important part of pathogen detection moving forward.

## Citing the Workflow

A manuscript describing these results and the experimental design that produced them is in preparation. We will update this section when a preprint and the eventual publication is available.
