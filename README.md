# Functional Annotation Pipeline

This pipeline is designed to perform comprehensive functional annotation of genomic data, specifically focusing on proteins. It integrates various bioinformatics tools to predict protein functions, domains, subcellular localization, and other important features.

## Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Input Files](#input-files)
- [Configuration](#configuration)
- [Usage](#usage)
- [Pipeline Workflow](#pipeline-workflow)
- [Output Files](#output-files)
- [Troubleshooting](#troubleshooting)
- [Acknowledgments](#acknowledgments)
- [License](#license)
- [Contact](#contact)

## Overview

The Functional Annotation Pipeline automates the process of annotating proteins from genomic data. By integrating tools like MMseqs2, InterProScan, HMMER, SignalP, TargetP, and MultiLoc2, it provides a streamlined workflow for:

- Filtering and processing GFF files
- Extracting protein and coding sequences
- Searching protein databases for homologs
- Predicting protein domains and functions
- Determining subcellular localization
- Identifying signal peptides and targeting signals
- Merging all annotations into a comprehensive report

## Features

- **Automated Workflow**: Utilizes Snakemake for workflow management.
- **Containerized Tools**: Employs Singularity containers for reproducibility.
- **Customizable**: Configurable via a YAML file.
- **Scalable**: Supports multithreading and parallel processing.
- **Comprehensive Output**: Generates a detailed TSV file with annotations.

## Prerequisites

- **Snakemake**: Workflow management system
- **Singularity**: Containerization platform

### Singularity Images

The necessary `.sif` and `.img` Singularity images are hosted on [Sylabs Cloud](https://cloud.sylabs.io/library/danifilho/functional_annotation_images). You can download them using the following commands:

```bash
singularity pull agat.sif library://danifilho/functional_annotation_images/agat:latest 
singularity pull gffread.sif library://danifilho/functional_annotation_images/gffread:latest 
singularity pull hmmer3.sif library://danifilho/functional_annotation_images/hmmer3:latest
singularity pull iprscan.sif library://danifilho/functional_annotation_images/iprscan:latest
singularity pull mmseqs2.sif library://danifilho/functional_annotation_images/mmseqs2:latest
singularity pull multiloc2_v3.img library://danifilho/functional_annotation_images/multiloc2_v3:latest
singularity pull samtools.sif library://danifilho/functional_annotation_images/samtools:latest
singularity pull sigtarp.sif library://danifilho/functional_annotation_images/sigtarp:latest
```

## Installation

**Ensure Singularity is installed**:

   - Install Singularity following the official [Singularity Installation Guide](https://sylabs.io/guides/3.0/user-guide/installation.html).

## Dependencies

### Software

- **Snakemake**: Workflow management system
- **Singularity**: Containerization platform

### Python Packages

- **biopython**

## Input Files

Place the following input files in the `inputs/` directory:

- **Genome FASTA File**: e.g., `genome.fasta`
- **GFF Annotation File**: e.g., `annotations.gff`

Place the custom scripts in the `scripts/` directory:

- `functional_merge.py`
- `multiloc_script.py`

## Configuration

Create a `config.yaml` file in the root directory with the following content:

```yaml
species_name: "your_species_name"
volume_name: "/absolute/path/to/your/working/directory"
fasta_file: "genome.fasta"
final_filtering_gff: "annotations.gff"
multiloc_script: "multiloc_script.py"
mmseqs_databases:
  - name: "database1"
    path: "/absolute/path/to/databases/database1"
  - name: "ncbi"
    path: "/absolute/path/to/ncbi_data_file.prot"
```

**Important**: Use absolute paths for `volume_name` and database paths.

## Usage

### Running the Pipeline

```bash
snakemake --cores <number_of_cores> --use-singularity
```

Replace `<number_of_cores>` with the number of CPU cores available.

### Dry Run

To perform a dry run and check the workflow:

```bash
snakemake -n
```

### Cleaning Up

To clean up all generated files:

```bash
snakemake clean
```

## Pipeline Workflow

The pipeline consists of the following steps:

1. **All**: Specifies the final output.
2. **AGAT**: Filters the GFF file to keep the longest isoform.
3. **GFFRead**: Extracts protein (AA) and coding sequences (CDS).
4. **MMseqs2**: Searches protein sequences against specified databases.
5. **Best Hits**: Extracts the best hit from MMseqs2 outputs.
6. **InterProScan**: Performs functional annotation.
7. **Protein List Generation**: Prepares inputs for MultiLoc2.
8. **MultiLoc2**: Predicts subcellular localization.
9. **HMMER**: Identifies protein domains using Pfam-A.
10. **SignalP and TargetP**: Predicts signal peptides and targeting signals.
11. **Generate TSV**: Merges all annotations into a TSV file.

## Output Files

All output files are placed in their respective directories:

- **AGAT Output**: `agat_outputs/`
- **GFFRead Output**: `gff_read_output/`
- **MMseqs2 Output**: `mmseqs_output/`
- **InterProScan Output**: `iprscan_output/`
- **MultiLoc2 Output**: `multiloc_outputs/`
- **HMMER Output**: `hmmer_outputs/`
- **SignalP and TargetP Output**: `sigtarp_outputs/`
- **Final Annotation File**: `functional_outputs/functional_annotation.tsv`

- Output Structure:
ID sequence: The unique identifier for each gene.

NCBI DESCRIPTION: Description of the gene from the NCBI database.

NCBI Subject; E-value; Bit score: NCBI annotation details, including the subject identifier, E-value, and bit score for the gene.

Other MMseqs database columns (dynamic): Each additional MMseqs database used in the analysis will have a set of columns named {Database} Subject; E-value; Bit score, containing the corresponding annotation details for that database.

SignalP Pos; Pr: SignalP data with position and prediction probability.

TargetP Prediction; noTP; SP; mTP; cTP; luTP; CS Position: TargetP data with localization prediction and position details.

MULTILOC: MultiLoc localization prediction.

IPRSCAN GO: Gene Ontology (GO) terms from IPRScan.

IPRSCAN IPR: InterPro (IPR) annotations from IPRScan.

Hmmer Pfam: Pfam domain annotations from HMMER.

Description
The TSV file is a comprehensive summary of gene annotations, dynamically incorporating columns for multiple MMseqs databases. Each row represents one gene, integrating data from several annotation sources (NCBI, SignalP, TargetP, MultiLoc, IPRScan, and HMMER), providing a unified view of functional characteristics across databases.

## Troubleshooting

- **Singularity Issues**: Ensure that paths in `config.yaml` are absolute and that Singularity has the necessary permissions.
- **File Not Found Errors**: Verify that all input files and databases are correctly placed and paths are accurate.
- **Resource Limitations**: Adjust `--cores` based on your system's capabilities.
- **Permission Errors**: Check file and directory permissions, especially when writing outputs.

## Acknowledgments

This pipeline utilizes open-source tools and databases. We acknowledge the developers and maintainers of these resources:

- **AGAT**
- **gffread**
- **MMseqs2**
- **InterProScan**
- **HMMER**
- **SignalP**
- **TargetP**
- **MultiLoc2**

## License

(LICENSE)

## Contact

For any questions or issues, please contact:

- **idk yet what to put here**
- **Email**: email@example.com
- **GitHub**: [username](https://github.com/username)

