#!/bin/bash
#
# Script Name: 5_Hisat2.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script edits raw fastq files and once those are ready, it aligns them to reference genome in
# also marks PCR duplicates in your alignment file.
#
# Run Information: This script is run manually by '5_Hisat2.sh'.
#
# Standard Output: The output is mapped files.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirenment: Prepare a reference genome before alignment, download fasta (*fa.gz) file that contains DNA sequence for S. pombe https://fungi.ensembl.org/Schizosaccharomyces_pombe/Info/Index
# once you download it, copy it to workdir/reference
#


cd ~/workdir/
source activate fastq

# Move the UMIs (Unique Molecular Identifier) to headers, here it moves 8 bases, to change it UMI1:<put number of bases to be moved here>

for file in *.fastq
do
	je clip F=$file 'RL=<UMI1:8><SAMPLE1:x>' 'OL=1:U1:S1' OUTPUT_DIR=edit #je clip .....
	gunzip edit/out_1.txt.gz
	mv edit/out_1.txt edit/c_$file
	gzip $file
done
ls -l edit/ 

# For each read remove the first nucleotide that comes from ligation step (-f 2) and the last one that usually has odd base distribution (-t 1).
cd edit/
for file in c_*.fastq
do
	fastx_trimmer -f 2 -i $file -o t_$file
	rm $file
	fastx_trimmer -t 1 -i t_$file -o tt_$file
	rm t_$file
done
ls -l

# Add index to the reference genome
cd ../reference/
gunzip *.gz
hisat2-build *.fa reference

# This part divides the file into several smaller files named "reference?.ht2". The files are ready to be aligned.

cd ../edit/
for file in *_R1.fastq # make the iteration of files "*_R1.fastq"
do 
	base=${file%_R1*} # delete the "_R1" and another part from the file name (remain "*")
	x=${base}_R1.fastq y=${base}_R2.fastq; # create two variables: "*_R1.fastq" and "*_R2.fastq"
	echo $x $y # show on the screen
	short=${base#*_c_} # It deletes all until "_c_" (including "_c_"), and save it in varianble "short". Example: "read_c_numberone" >> "numberone"
	echo $short
	hisat2 -q -x ../reference/reference -1 $x -2 $y --rna-strandness FR -S ../mapped/${short}.sam --summary-file  ../mapped/hisat2/${short}.txt --new-summary -p 16
	rm $x # remove "*_R1.fastq" and "*_R2.fastq"
	rm $y
done

# -q indicates that input file is in .fastq format,
# -x ../prefix indicates which files will be used as the reference,
# -S is used to define the output file
#
# Introducing "base" is necessary to pair reads.
# "%_" shorten name by cutting off everything that stands after %
#
# I also introduced "short" to make the .sam file look like <genename.sam>
# "#*_c_" shorten name by cutting off everything that stands pre "_c_" including this motive


# Sort .sam files with samtools
# Mark PCR duplicates with "je" and check out output file with statistics.

cd ../mapped
for file in *.sam
do
	samtools sort $file -o s_$file -@ 8
	rm $file
done

for file in s_*.sam
do
	base=${file#s_} # name of the file without "s_" at the beginning
	red='\033[0;31m' # bash red color
	white='\033[0m'
	printf "${red}$base${white}"
	short=${base%.sam} #name of the file without ".sam" at the end
	je markdupes I=$file O=m_$file M=markdupes/dups_report_${short}.txt MM=0 
	rm $file 
done
