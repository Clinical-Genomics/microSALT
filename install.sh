#!/usr/bin/env bash

set -e
set -x
shopt -s nullglob

#Suggests provided branch. Else suggests master
default_branch=${1-master}

echo "Welcome to the microSALT installation script. Q to exit"
while true; do
    echo "Would you like a 'release' or 'source' (development) environment ['release']?"
    read input
    if [[ $input = "q" ]] || [[ $input = "Q" ]]; then
        exit 0
    elif [[ $input = "y" ]] || [[ $input = "yes" ]] || [[ $input = "" ]]; then
        type="release"
        break
    elif [[ $input == "source" ]] || [[ $input == "release" ]]; then
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
    rm -rf microSALT
fi

git clone https://github.com/Clinical-Genomics/microSALT
cd microSALT && git checkout $branch

if [[ $type == "release" ]]; then
    uv sync
elif [[ $type == "source" ]]; then
    uv sync --group dev
fi

echo "Installation Complete! Activate the environment with: source microSALT/.venv/bin/activate"
echo "Configuration requires manual set-up as described in README.md"
