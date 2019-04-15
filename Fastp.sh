#!/bin/bash
#
# Script Name: Fastp.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script merges FastQ files (the reads R1 with R2, respectively), in case of their high 
# similarity or low quality of at least one of them.
#
# Run Information: This script is run manually by 'bash Fastp.sh'.
#
# Standard Output: The output is a file merging two reads into one. 
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirement: Install fastp according to instruction at https://github.com/OpenGene/fastp, for example with bioconda:
#'conda install -c bioconda fastp'
#


for files in *R1.fastq
do
base=${files%_R1*}
echo $base
x=${base}_R1.fastq
y=${base}_R2.fastq
echo $x $y
#short=${base#_c_}
#echo $short
fastp -i $x -o ${short}.fastq [-I $y] [-m][--discard_unmerged][--overlap_len_require 10][-j fastp/${short}.json][-h fastp/${short}.html][-A]

# -A, --disable_adapter_trimming, adapter trimming is enabled by default.
