#!/bin/bash
#
# Script Name: Cutadapt.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script can be implemented for the reads from Illumina that are longer than library preps. 
# It cuts off adapters sequences and everything that follows them.
#
# Run Information: This script is run manually by 'bash Cutadapt.sh'.
#
# Standard Output: The output is a fastq file with removed adaptors.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirement: Install cutadapt with other requirements or using "pip3 install --user --upgrade cutadapt"; 
# Download BBMap_38.42.tar.gz file to your home directory (https://sourceforge.net/projects/bbmap/files/) (or the other version).
#

cd ~/workdir/
source activate fastq

for files in *R2.fastq
do
base={file%.fastq}
cutadapt -a GATCGTCGGACTGTAGAACTCTGAAC -o edit/t_$files $files -m 25 --info-file edit/report_$base.txt
done 

# The cutadapt command recognizes Illumina TruSeq Small RNA 5' adapter (given sequence is reverse complement).
#
# -o is output; input follows it without any flag,
# -m option saves only those reads that remain longer than 25 nucleotides after cutting,
# --info-file is used to save report file.
#
# You can also use -j option that allows parallel runs by using multiple cores (-j N), N stands for number of cores.

for files in *R1.fastq
do
base={file%.fastq}
cutadapt -a TGGAATTCTCGGGTGCCAAGG -o edit/t_$files $files -m 25 --info-file edit/report_$base.txt
done 


# In our case, *R2.fastq reads had 3 nucleotides added at the beginning of the sequence (we don't know the reason).
# You can remove them with other cutadapt command or with fastx_trimmer.

cd edit/
for files in *R2.fastq
do
cutadapt -u 3 -o t$files $files
done

# If you want to remove nucelotides standing at the end of the sequence add a minus sign (e.g. -u -3).
#
# Before aligning make sure before aligning the reads that each read in R1 file has its pair in R2 file.
# You can use bbmap package which cant be downloaded by conda :(
#cd
#get BBMap_38.42.tar.gz file to your home directory (https://sourceforge.net/projects/bbmap/files/) (or other version that is more up to date)
#tar -xvzf BBMap_(version).tar.gz
#now use repair.sh script to get rid off unpaired reads, at the same time u can rename files so that u can put them to 3_hisat2 script

cd ~/workdir/edit
tar -xvzf BBMap_38_42.tar.gz
for file in *_R1.fastq; do  base=${file%_R1*}; echo $base; x=${base}_R1.fastq y=t${base}_R2.fastq; echo $x $y; short=${base#*t_}; echo $short; ../../bbmap/repair.sh in1=$x in2=$y out1=$tt_c_{short}_R1.fastq out2=$tt_c_{short}_R2.fastq outs=bbtools/singletons_${short}.fastq repair ; done
