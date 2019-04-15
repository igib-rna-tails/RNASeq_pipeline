#!/bin/bash
#
# Script Name: 7_IGV.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script prepares your mapped files for analysis in IGV software.
#
# Run Information: This script is run manually by '7_IGV.sh'.
#
# Standard Output: The output is the set of mapped files for analysis in IGV software.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirement: Installation of IGV, virtual environment, and files "m_s_*.sam" created at the previous steps.
#

cd ~/workdir/mapped
source activate fastq

# Convert your .sam files to .bam
for files in m_s_*.sam
do
	base=${files%.sam}
	samtools view -Sb $files -o ${base}.bam
done

# Add index to the BAM files
#
#-b generates .bai file
#
for files in *.bam
do
	samtools index -b $files
done

# Start analyse your mapping in IGV
source activate igv
igv
