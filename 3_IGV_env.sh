#!/bin/bash
#
# Script Name: 3_IGV_env.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script creates the environment for IGV, which can't be installed using fastq environment.
#
# Run Information: This script is run manually by 'bash 3_IGV_env.sh'.
#
# Standard Output: The output is a new environment for IGV.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirenment: Download of Miniconda for your operating system from https://conda.io/miniconda.html
#


# Create an another environment for IGV (it can't be installed to fastq env)
conda create --name igv
conda activate igv
conda install -c bioconda igv
conda deactivate igv
