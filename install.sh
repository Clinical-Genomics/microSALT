#!/usr/bin/env bash

set -e
shopt -s nullglob

#Suggests provided branch. Else suggests master
default_branch=${1-master}
default_name=${2-microSALT}

echo "Welcome to the microSALT installation script. Q to exit"
while true; do
    echo "Name your microSALT environment ['microSALT']:"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input = "y" ]] || [[ $input = "yes" ]] || [[ $input = "" ]]; then
        cname=$default_name
        break
    else
        cname=$input
        break
    fi
done
while true; do
    echo "Would you like a 'release' or 'source' (development) environment ['release']?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        break
    elif [[ $input = "y" ]] || [[ $input = "yes" ]] || [[ $input = "" ]]; then
        type="release"
        break
    elif [[ $input == "source" ]]  || [[ $input == "release" ]]; then
        type=$input
        break
    fi
done

validbranch=false
while true; do
    echo "Name the branch to install ['$default_branch']:"
    while ! $validbranch; do
        read input
        if [[ $input = "y" ]] || [[ $input = "yes" ]] || [[ $input = "" ]]; then
            branch=$default_branch
        else
            branch=$input
        fi
        curl https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/LICENSE | tac | tac | grep -q 'License' && validbranch=true||echo "Invalid branch name"
    done
    break
done
echo "Thank you, setting up environment $cname!"

#Unload environment
conda info | tac | tac | grep -q $cname && source deactivate || :
#Remove environment if already present
conda remove -y -n $cname --all || :

conda create -y -n $cname python=3.6
conda activate $cname
echo "add defaults"
conda config --add channels defaults
echo "add bioconda"
conda config --add channels bioconda
echo "add conda-forge"
conda config --add channels conda-forge
echo "add blast"
conda install blast=2.12.0
echo "add bwa"
conda install bwa=0.7.17
echo "add picard"
conda install picard=2.20.3
echo "add pigz"
conda install pigz=2.4
echo "add quast"
conda install quast=5.0.2
echo "add samtools"
conda install samtools=1.13
echo "add spades"
conda install spades=3.13.1
echo "add trimmomatic"
conda install trimmomatic=0.39
echo "add r-base"
conda install r-base=4.1.1

if [[ $type == "release" ]]; then
    pip install -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/requirements.txt -r https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/requirements-dev.txt 
    pip install -U git+https://github.com/Clinical-Genomics/microSALT@$branch
elif [[ $type == "source" ]]; then
  HERE=$PWD
  if [ -d ${HERE}/microSALT ]; then
    rm -rf microSALT
  fi
  git clone https://github.com/Clinical-Genomics/microSALT
  cd microSALT && git checkout $branch
  pip install -r requirements.txt -r requirements-dev.txt && pip install -e . && cd ${HERE}
  echo "Source installed under ${HERE}/microSALT" 
fi 
echo "Installation Complete!"
while true; do
    echo "Configuration requires manual set-up as described in README.md ['yes']:"
    read input
    if [[ $input = "y" ]] || [[ $input = "yes" ]]; then
        break
    fi
done
