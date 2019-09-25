#!/bin/bash
#
# Script Name: 6_RSeqC.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script allows to check theoritical distances between paired reads using RSeqC.
#
# Run Information: This script is run manually by '6_RSeqC.sh'.
#
# Standard Output: The output is *gff3.bed file.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirement: Download annotated reference genome with the extension .gff3 (or similar). If it differs change code # below. 
# Move the file to ~/workdir/reference2/.
#


source activate fastq

current=`pwd` #current working directory
cd ~/workdir/reference2
gff2bed < *.gff3 > reference.gff3.bed #If you prepare annotated reference genome with differ extension (not .ggf3), please change this part of code.

# RSeqC will create 4 files, one of which is run with R
#
#-i is input, -r is reference, -o is output

cd ../mapped
for file in m_s_*.sam
do
	base=${file#*m_s_}
	inner_distance.py -i $file -r ../reference2/reference.gff3.bed -o rseqc/report_$base
done

cd ${current}
