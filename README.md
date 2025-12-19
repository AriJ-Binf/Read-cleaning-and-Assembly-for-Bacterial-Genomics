# Read Cleaning and Assembly for Bacterial Genomics
A flexible, end-to-end read trimming, QC, and genome assembly pipeline for paired-end Illumina data.
Supports Trimmomatic or fastp for read trimming and performs SPAdes assembly for each sample automatically.

Designed for bacterial genomics workflows using SRA-style paired FASTQ files.

Overview

This script performs the following steps for each paired-end sample in an input directory:

Quality Control

FastQC (Trimmomatic mode only)

fastp HTML/JSON reports (fastp mode)

Read Trimming

Adapter trimming (optional)

Quality trimming (sliding window, leading/trailing)

Minimum read length filtering

Singleton read handling

Genome Assembly

SPAdes (paired reads + optional singletons)

Assembly logs and contigs output per sample

Each sample is processed independently and written to its own output directory.

Input Requirements

Paired-end FASTQ files named:

SAMPLE_1.fastq.gz
SAMPLE_2.fastq.gz


All input files must be located in a single input directory

Gzipped FASTQ files are required

Directory Structure (Output)

For each sample (SAMPLE), the following structure is created:

SAMPLE/
├── QC/
│   └── FastQC reports (Trimmomatic mode only)
├── trim/
│   ├── r1.paired.fq.gz
│   ├── r2.paired.fq.gz
│   ├── singletons.fq.gz
│   ├── trimming logs
│   └── fastp.html / fastp.json (fastp mode)
└── asm/
    └── spades/
        ├── contigs.fasta
        ├── spades.stdout.txt
        └── spades.stderr.txt

Dependencies

Make sure the following tools are available in your PATH:

bash (≥ 4)

FastQC

Trimmomatic or fastp

SPAdes

gzip

This script is compatible with conda / mamba environments and macOS or Linux systems.

Usage
bash Trim_and_Assembly.sh [OPTIONS] <raw_data_folder>

Example (Trimmomatic)
bash Trim_and_Assembly.sh \
  --trimmer trimmomatic \
  --threads 16 \
  --phred 33 \
  --slidingwindow 5:30 \
  --minlen 50 \
  raw_fastq/

Example (fastp)
bash Trim_and_Assembly.sh \
  --trimmer fastp \
  --threads 16 \
  --detect-adapters \
  raw_fastq/
