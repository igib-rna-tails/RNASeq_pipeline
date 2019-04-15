#!/bin/bash
#
# Script Name: 4_Fastqc.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The script activates fastq env, makesgunzipes all of the *fastq.gz files, and make their quality control with FastQC.
#
# Run Information: This script should run manually by 'bash 4_Fastqc.sh'.
#
# Standard Output: The output is the zipped fastqc in ~/workdir/fastqc/.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirement: Make a copy of the *.fastq or *.fastq.gz files in ~/workdir/.
#


cd ~/workdir/

source activate fastq

gunzip *.gz

# Make quality control with fastqc
#
# -o stands for output
#
for fastqfiles in *.fastq

do
	fastqc $fastqfiles -o fastqc/
done
