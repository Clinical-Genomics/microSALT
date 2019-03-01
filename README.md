<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

__Microbial Sequence Analysis and Loci-based Typing pipeline__

_The microbial sequence analysis and loci-based typing pipeline (microSALT) is used to analyse microbial samples.
It produces a quality control of the sample, determines a sample's organism specific sequence type, and its resistance pattern. microSALT also provides a database storage solution and report generation of these results.

microSALT uses a combination of python, sqLite and flask. Python is used for the majority of functionality, the database is handled through sqLite and the front-end is handled through flask. All analysis activity by microSALT requires a SLURM cluster._

## Quick installation
```
conda config --add channels bioconda
conda create -n microSALT python=3.6
source activate microSALT
conda install -c bioconda blast=2.5.0=h3727419_3 spades=3.12.0=py36_0 \
trimmomatic=0.38=1 bwa==0.7.15=1 samtools=1.6=0 picard=2.18.26=0
git clone https://github.com/Clinical-Genomics/microSALT.git
cd microSALT && pip install -r requirements.txt && pip install.
```

Then continue with __Configuration__

## Quick installation
* `git clone https://github.com/Clinical-Genomics/microSALT.git`
* `cd microSALT && bash install.sh`
* Perform all steps under section  __Configuration__

## Configuration
Copy the configuration file to microSALTs hidden home directory, _or_ copy the configuration file anywhere and direct the envvar MICROSALT_CONFIG to it. See examples: 

`cp configExample.json $HOME/.microSALT/config.json`

_or_
```
cp configExample.json /MY/FAV/FOLDER/config.json
export MICROSALT_CONFIG=/MY/FAV/FOLDER/config.json
```

__Then edit the fields to match your environment__.

### LIMS Configuration
Create `$HOME/.genologicsrc` with the following formatting:
```
[genologics]
BASEURI=https://yourlims.corporation.se/
USERNAME=your_username
PASSWORD=your_password
[logging]
MAIN_LOG=/tmp/lims.log
```

## Usage
* Use the `analyse` function to start sbatch job(s), producing output to `folders['results']`. Afterwards the parsed results  are uploaded to the SQL back-end and produce reports (HTML).
* Various functionality, including adding manually new reference organisms and re-generating reports; are stored under the `utils` commands.

## Databases
### MLST Definitions
microSALT will automatically download & use the MLST definitions for any organism on pubMLST (https://pubmlst.org/databases/).
Other definitions may be used, as long as they retain the same format. 

### Resistance genes
microSALT will automatically download & use the resistance genes of resFinder (https://cge.cbs.dtu.dk/services/data.php).
Any definitions will work, as long as they retain the same formatting.

## Requirements
### Hardware
* A slurm enabled HPC
* A (clarity) LIMS server
* A sqLite service

### Software
* Conda
* Python 3.6
