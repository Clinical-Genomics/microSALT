#!/bin/bash

set -e
shopt -s nullglob

echo "Welcome to the microSALT installation script. Q to exit"
echo "What name would you like to give your conda environment?"
while true; do
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]] 
        then break 
    else 
        cname=input
    fi
done
echo "Would you like a 'dev' or 'prod' environment?"
while true; do
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]
        then break
    else
        type=input
    fi
done
echo "Thank you, executing!"

conda create -n $cname python=3.6
source activate $cname
conda config --add channels bioconda
conda install -c bioconda blast=2.5.0=hc0b0e79_3 spades=3.12.0=py36_0 trimmomatic=0.38=1
git clone https://github.com/Clinical-Genomics/microSALT.git
if $type = 'prod'
    cd microSALT && pip install -r requirements.txt && pip install .
elif $type = 'dev'
    cd microSALT && pip install -r requirements.txt && pip install -e .

echo "Installation complete. Remember to set-up a configuration as mentioned in the README.md"
