#!/bin/bash
#
# Script Name: 1_Conda.sh
#
# Authors: Maciej Grochowski, Michal Malecki, Lidia Lipinska (Documentation)
# Date: 28.03.2019
#
# Description: The following script downloads miniconda.
#
# Run Information: This script is run manually by 'bash 1_Conda.sh' (conda installation), and next commands manually in bash.
#
# Standard Output: The output is an installation of conda, necessary packages, and creation of new environment for IGV.
#
# Error Log: No existing file to storage the errors associated with this script.
#
# Requirenment: Download Miniconda for your operating system from https://conda.io/miniconda.html
#


# Please download Miniconda and run following command to install it.
bash ~/Pobrane/Miniconda3-latest-Linux-x86_64.sh

# If necessary, restart terminal.
exec bash
read -p "Press [Enter] key to start backup..."
conda list


