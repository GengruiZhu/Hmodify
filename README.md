# Hmodify
# Chromosome Manipulation Pipeline

A bioinformatics pipeline for manipulating chromosome structures based on AGP and FASTA files, designed for Hi-C scaffolding applications.

## Description

This tool performs complex chromosome operations including:
- Extracting chromosome segments from AGP files
- Copying segments between chromosomes
- Inserting segments at specific locations
- Generating modified FASTA files reflecting the new chromosome structures

## Features

- **Flexible Configuration**: Define operations via a simple config file
- **Logging**: Detailed logging for troubleshooting
- **Validation**: Input validation and error checking
- **Modular Design**: Each operation is performed in discrete steps
- **Reference Support**: Optional reference chromosome integration

## Requirements

- Python 3.6+
- seqkit (for FASTA extraction)
- Standard UNIX tools (awk, grep, cat, etc.)

## Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/GengruiZhu/Hmodify.git
   cd Hmodify
