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
    echo "Would you like a 'dev' or 'prod' environment?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input != "dev" ]] && [[ $input != "prod" ]]; then
        :
    else
        type=$input
        break
    fi
done
if [[ $type = "prod" ]]; then
    rname=P_$cname
elif [[ $type = "dev" ]]; then
    rname=D_$cname
fi
echo "Thank you, setting up environment $rname!"

conda create -y -n $rname python=3.6
#source deactivate
source activate $rname
conda config --add channels bioconda
conda install -y -c bioconda blast=2.5.0=hc0b0e79_3 spades=3.12.0=py36_0 trimmomatic=0.38=1 samtools=1.6=0 picard=2.18.26=0
if [[ $type = "prod" ]]; then
    pip install -r requirements.txt && pip install .
elif [[ $type = "dev" ]]; then
    pip install -r requirements.txt && pip install -e .
fi
echo "Installation Complete!"
while true; do
    echo "Write 'ok' to promise that you'll set-up the configuration as mentioned in README.md"
    read input
    if [[ $input = "ok" ]] || [[ $input = "OK" ]]; then
        break
    fi
done
