# Pipeline for analysis and visualization of data from RNA Sequencing
***
Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)

University of Warsaw, Institute of Genetics and Biotechnology
***
This script contains a procedure to analyze and visualize data from RNA Sequencing.
***
__The pipeline includes the following steps:__

    1. Installation of requirements.
    2. Creating the directories (folders) to save input and output of the analysis.
    3. Installation of Miniconda.
    4. Preparation of environment for conda.
    5. Preparation of environment for IGV.
    6. Quality control using FastQC.
    7. Alignment using HISAT2.
    8. The quality control of RNA-seq using RSeQC.
    9. High-performance visualization of genomic datasets using IGV.
    10. Parsing with HTSeq.
    11. Checking the results carefully.
    12. If needed - usage of Fastp or Cutadapt to minimalize the bias.
*** 
    
### Requirements:
    * FastQC
    * MultiQC
    * BBmap
    * IGV
    * HISAT2
    * RSeQC
    * HTSeq
    * fastp
    * cutadapt
    
***
## Installation of requirements
#### A) Installation of basic requirements

You can install requirements using __*pip install*__ command. We recommend installing the requirements in a virtual environment to avoid inconsistencies in the local environment.

> pip install -r requirements.txt

#### B) Installation of FastQC
Install FastQC according to [this instruction](https://raw.githubusercontent.com/s-andrews/FastQC/master/INSTALL.txt).

#### C) Installation of BBmap
Install BBmap according to [this instruction](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/installation-guide/).

#### D) Installation of HISAT2
Install HISAT2 with conda and uptade it
> conda install hisat2
> conda update hisat2

#### E) Installation of HTSeq
Install HTSeq according to [this instruction](https://htseq.readthedocs.io/en/release_0.11.1/install.html).
Installation with "sudo apt-get". At first, check your version of Python and install HTSeq for your version.
> python --version

##### Choose your Python version
> sudo apt-get install build-essential python3.7-dev python-numpy python-matplotlib python-pysam python-htseq
#### F) Installation of Fastp with Bioconda
> conda install -c bioconda fastp 
