#!/bin/bash
#
# Script Name: 8_HtSeq.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script prepares SAM files and creates a file, that shows how many reads are mapped 
# to any particular feature annotated in reference genome.
#
# Run Information: This script is run manually by '8_HtSeq.sh'.
#
# Standard Output: The output is a file, that shows how many reads are mapped to any particular feature annotated in reference genome.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirenment: Prepared fastq environment, installed requirements, and "m_s_*.sam" created in the previous steps.
#

# Prepare the SAM file to make HtSeq analysis

cd ~/workdir/mapped
source activate fastq

# Sort the reads by names
for file in m_s_*.sam
do
	base=${file#m_s_}
	samtools sort -n $file -o sbn_$base
	rm $file
done


# By default feature=exon (use -t <featurename> option to change it).
# Script searches for feature name in third column of .gff3 file,then jumps to 8th column and by deafult looks for gene_id,
# which might be not present there. In such a case use -i <8thcolumnstartingword> option
#
# HINT: REMEMBER to change code if your reference has different extension than .gff3.
# 
# Depending on the method of library preparation Read1 might be equal to feature sequence or complementar to it.
# By default its treated as equal, # here we change it with -s <reverse> option.


for files in sbn_*.sam
do
	base=${files%.sam}
	short=${base#sbn_}
	red='\033[0;31m'
	white='\033[0m'
	printf "${red}$short${white}" 
	htseq-count -t gene -i ID -s reverse $files ../reference2/*.gff3 > htseq/counttable_gene_${short}.txt
done

# The command below prints the sum of every score located in second column of the .txt file, except last five lanes.

cd htseq
for file in *.txt
do
	red='\033[0;31m'
	white='\033[0m'
	printf "${red}$file${white}" 
	head -n -5 $file | cut -f 2 | awk '{total += $0} END{print "sum="total}'
done
