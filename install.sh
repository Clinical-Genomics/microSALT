#!/usr/bin/env bash

set -e
shopt -s nullglob

echo "Welcome to the microSALT installation script. Q to exit"
while true; do
    echo "What name would you like to give your microSALT environment (no prefix)?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    else
        cname=$input
        break
    fi
done
while true; do
    echo "Prod = Latest release. Stage = Any branch. Dev = Source-modifiable installation"
    echo "Would you like a 'prod', 'stage' or 'dev' environment?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input == "prod" ]]; then
        break
    elif [[ $input == "dev" ]] || [[ $input == "stage" ]]; then
        type=$input
        validbranch=false
        while ! $validbranch; do
            echo "Name a branch to install (or just 'master')":
            read input
            branch=$input
            curl https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/LICENSE | grep -q 'License' && validbranch=true||echo "Invalid branch name"
        done
        break
    fi
done

if [[ $type = "prod" ]]; then
    rname=P_$cname
elif [[ $type = "dev" ]]; then
    rname=D_$cname
elif [[ $type = "stage" ]]; then
    rname=S_$cname
fi
echo "Thank you, setting up environment $rname!"

conda info| grep -q $rname && source deactivate || ""

#Accepts that environment doesnt exist
conda remove -y -n $rname --all || ""
conda create -y -n $rname python=3.6
#source deactivate
source activate $rname
conda config --add channels bioconda
conda install -y -c bioconda blast=2.9.0 bwa=0.7.17 picard=2.20.3 quast=5.0.2 samtools=1.9 spades=3.13.1 trimmomatic=0.39
if [[ $type = "prod" ]]; then
  pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/requirements.txt 
  pip install -U git+https://github.com/Clinical-Genomics/microSALT 
elif [[ $type = "stage" ]]; then
  pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/requirements.txt &&
  pip install -U git+https://github.com/Clinical-Genomics/microSALT@$branch
elif [[ $type = "dev" ]]; then
  pip install -r requirements.txt && pip install . && python setup.py develop
fi 
echo "Installation Complete!"
while true; do
    echo "Write 'ok' to promise that you'll set-up the configuration as mentioned in README.md"
    read input
    if [[ $input = "ok" ]] || [[ $input = "OK" ]]; then
        break
    fi
done
