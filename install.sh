#!/usr/bin/env bash

set -e
set -x
shopt -s nullglob

#Suggests provided branch. Else suggests master
default_branch=${1-master}

echo "Welcome to the microSALT installation script. Q to exit"
while true; do
    echo "Would you like a 'production' or 'stage' (development) environment ['production']?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        exit 0
    elif [[ $input = "y" ]] || [[ $input = "yes" ]] || [[ $input = "" ]]; then
        type="production"
        break
    elif [[ $input == "stage" ]] || [[ $input == "production" ]]; then
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
        curl https://raw.githubusercontent.com/Clinical-Genomics/microSALT/$branch/LICENSE | tac | tac | grep -q 'License' && validbranch=true || echo "Invalid branch name"
    done
    break
done
echo "Thank you, installing branch $branch!"

if [ -d microSALT ]; then
    while true; do
        echo "Directory microSALT already exists. Do you want to delete it and continue? [y/N]"
        read input
        if [[ $input = "y" ]] || [[ $input = "yes" ]]; then
            rm -rf microSALT
            break
        elif [[ $input = "n" ]] || [[ $input = "no" ]] || [[ $input = "" ]]; then
            echo "Exiting installation"
            exit 0
        fi
    done
fi

git clone https://github.com/Clinical-Genomics/microSALT
cd microSALT && git checkout $branch

if [[ $type == "production" ]]; then
    uv sync
elif [[ $type == "stage" ]]; then
    uv sync --group dev
fi

echo "Installation Complete! Activate the environment with: source microSALT/.venv/bin/activate"
echo "Configuration requires manual set-up as described in README.md"
