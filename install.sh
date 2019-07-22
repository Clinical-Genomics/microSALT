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

conda create -y --force -n $rname python=3.6
#source deactivate
source activate $rname
conda config --add channels bioconda
conda install -y -c bioconda blast=2.9.0=pl526h979a64d_3 bwa=0.7.17=h84994c4_5 picard=2.20.3=0 \
quast=5.0.2=py36pl526ha92aebf_0 samtools=1.9=h8571acd_11
pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/requirements.txt 
pip install -U git+https://github.com/Clinical-Genomics/microSALT 
echo "Installation Complete!"
while true; do
    echo "Write 'ok' to promise that you'll set-up the configuration as mentioned in README.md"
    read input
    if [[ $input = "ok" ]] || [[ $input = "OK" ]]; then
        break
    fi
done
