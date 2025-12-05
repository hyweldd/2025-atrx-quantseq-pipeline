# 2025-atrx-quantseq-pipeline
Code used to create interferon genes heatmap from QuantSeq data

The correct software versions for this code can be found in *bin/conda_create_command.sh*. 
This script can be used to create a conda environment containing the appropriate sofware. 
The software has been tested on Linux and Mac OSX. No non-standard hardware is required.

## Installation guide 

To install the software dependencies for nextflow and snakemake:
1. Run the *bin/conda_create_command.sh* using bash
2. Activate the conda environment

To install the software depencies for the custom R scripts:
1. Install R and *renv*
2. Initialise the *renv* environment using the *DESCRIPTION* file in the *Rproject* directory

The alignment pipeline uses igenomes for the reference data. The required files are specified in params.yaml.

## Usage

1. Apply the *bin/trim_all_polyA_plus_downstream_seq.py* script on the source *fastq.gz* files
2. Run the *bin/run_nfcore.sh* script to generate gene level count data 
3. Run *snakemake* using the snakefile *bin/Snakefile* to calculate differential expression 

