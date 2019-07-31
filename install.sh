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
    echo "Would you like a 'dev', 'stage' or 'prod' environment?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input == "dev" ]] || [[ $input == "prod" ]] || [[ $input == "stage" ]]; then
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
elif [[ $type = "dev" ]] || [[ $type = "stage" ]]; then
  pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/requirements.txt &&
  pip install -U git+https://github.com/Clinical-Genomics/microSALT@$branch
fi 
echo "Installation Complete!"
while true; do
    echo "Write 'ok' to promise that you'll set-up the configuration as mentioned in README.md"
    read input
    if [[ $input = "ok" ]] || [[ $input = "OK" ]]; then
        break
    fi
done
