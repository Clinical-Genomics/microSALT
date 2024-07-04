![Build status](https://github.com/clinical-genomics/microsalt/actions/workflows/run_tests.yml/badge.svg)
[![Coverage Status](https://coveralls.io/repos/github/Clinical-Genomics/microSALT/badge.svg?branch=master)](https://coveralls.io/github/Clinical-Genomics/microSALT?branch=master)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4026043.svg)](https://doi.org/10.5281/zenodo.4026043)

<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

__Microbial Sequence Analysis and Loci-based Typing pipeline__

_The microbial sequence analysis and loci-based typing pipeline (microSALT) is used to analyse microbial samples.
It produces a quality control of the sample, determines a sample's organism specific sequence type, and its resistance pattern. microSALT also provides a database storage solution and report generation of these results._

_microSALT uses a combination of python, sqLite and flask. Python is used for the majority of functionality, the database is handled through sqLite and the front-end is handled through flask. All analysis activity by microSALT requires a SLURM cluster._

## Quick installation
* `yes | bash <(curl https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/install.sh)`
* `cp configExample.json $HOME/.microSALT/config.json`
* `vim $HOME/.microSALT/config.json`

## Configuration
Copy the configuration file to microSALTs hidden home directory, _or_ copy the configuration file anywhere and direct the envvar MICROSALT_CONFIG to it. See example: 

`cp configExample.json $HOME/.microSALT/config.json`

_or_
```
cp configExample.json /MY/FAV/FOLDER/config.json
export MICROSALT_CONFIG=/MY/FAV/FOLDER/config.json
```

__Then edit the fields to match your environment__.

## Usage
* `microSALT analyse` contains functions to start sbatch job(s) & produce output to `folders['results']`. Afterwards the parsed results  are uploaded to the SQL back-end and produce reports (HTML), which are then automatically e-mailed to the user.
* `microSALT utils` contains various functionality, including generating the sample description json, manually adding new reference organisms and re-generating reports.

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


### Software
* Conda ( https://www.anaconda.com/distribution/ )
* Python 3.6
* sqLite Service ( https://www.sqlite.org/download.html )

## Contributing to this repo
This repository follows the Github flow approach to adding updates.
For more information, see https://guides.github.com/introduction/flow/

## Credits
* Isak Sylvin - Lead developer
* Emma Sernstad - Accreditation ready reports
* Tanja Normark - Various issues
* Maya Brandi - Various issues
