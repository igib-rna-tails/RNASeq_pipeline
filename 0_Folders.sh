#!/bin/bash
#
# Script Name: 0_Folders.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script creates all folders that are necessary for other scripts to work properly.
#
# Run Information: This script is run manually by 'bash 0_Folders.sh'.
#
# Standard Output: The output are the folders and sections to make following RNA-Seq Analysis.
#
# Error Log: No existing file to storage the errors associated with this script.
#

cd ~
mkdir workdir
cd workdir
mkdir fastqc
mkdir edit
mkdir mapped
mkdir reference
mkdir reference2
cd mapped/
mkdir hisat2
mkdir htseq
mkdir rseqc
mkdir markdupes
cd ~
