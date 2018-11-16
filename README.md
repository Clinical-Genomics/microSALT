<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

microbial Sequence Analysis and Loci-based Typing pipeline

The microbial sequence analysis and loci-based typing pipeline (microSALT) is used to determine a microbial sample's organism specific sequence type and resistance profile. This is in turn defined from a set of six to eight organism specific allele types, and a huge set of resistance genes. microSALT also provides a database storage solution and pdf generation of these results.

## Quick installation
### Conda dependency resolution
```
conda config --add channels bioconda
conda create -n microSALT python=3.6
source activate microSALT
conda install -c bioconda blast=2.5.0=h3727419_3 spades=3.12.0=py36_0 \
trimmomatic=0.38=1 bwa=0.7.17=ha92aebf_3 samtools=1.6=0
git clone https://github.com/Clinical-Genomics/microSALT.git
cd microSALT && pip install -r requirements.txt && pip install.
```
* Perform all steps under section  __Configuration__

## Configuration
Copy the configuration file `configExample.json` to `$HOME/.microSALT/config.json`
_or_
Copy the configuration file `configExample.json` to `/MY/FAV/FOLDER/config.json` 
and write `export MICROSALT_CONFIG=/MY/FAV/FOLDER/config.json`

Edit the fields to match your environment.

### LIMS Configuration
Create `$HOME/.genologicsrc` with the following formatting:
```
[genologics]
BASEURI=https://yourlims.corporation.se:8443
USERNAME=your_username
PASSWORD=your_password
[logging]
MAIN_LOG=/tmp/lims.log
```

## Usage
* Use the `start` function to start sbatch job(s), producing output to `folders['results']`. Afterwards the parsed results  are uploaded to the SQL back-end and produce reports (HTML).
* Various functionality, including adding manually new reference organisms and re-generating reports; are stored under the `util` command.

## Databases
### MLST Definitions
microSALT is able to neatly download the MLST definitions for any organism on pubMLST (https://pubmlst.org/databases/), and will attempt to automatically download these.
Other definitions may be used, as long as they retain the same format. 

### Resistance genes
microSALT relies on the resistance genes of resFinder (https://cge.cbs.dtu.dk/services/data.php), and will attempt to automatically download these.
Any definitions will work, as long as they retain the same formatting.

## Requirements
### Hardware
* A slurm enabled HPC
* A (clarity) LIMS server
* A sqLite service

### Software
* Conda
* Python 3.6
