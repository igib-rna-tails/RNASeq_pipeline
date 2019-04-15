#!/bin/bash
#
# Script Name: 2_Conda_env.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script creates the environment for conda packages and installs them.
#
# Run Information: This script is run manually by 'bash 2_Conda_env.sh'.
#
# Standard Output: The output is an installation of conda, necessary packages.
#
# Error Log: No existing file to storage the errors associated with this script.
#


# Create the environment for conda packages and install them in fastq env
conda create --name fastq
source activate fastq

conda install -c bioconda seqkit fastqc je-suite HISAT2 samtools bedops rseqc htseq cutadapt
#for packages in seqkit fastqc je-suite HISAT2 samtools bedops rseqc htseq cutadapt
#do
#	conda install -c bioconda $packages
#done

conda install -c biobuilds fastx-toolkit

source deactivate fastq

