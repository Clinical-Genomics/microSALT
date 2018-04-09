<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

microbial Sequence Analysis and Loci-based Typing pipeline

The microbial sequence analysis and loci-based typing pipeline (microSALT) is used to determine a microbial samples organism specific sequence type. This is in turn defined from a set of six to eight organism specific allele types. microSALT also provides a database storage solution and pdf generation of these results.

## Requirements
### Hardware
* A slurm enabled HPC
* A mySQL server
* A (clarity) LIMS server

### Software
* Conda
* Python 3.6

## Quick installation
### Conda dependency resolution
* `git clone https://github.com/Clinical-Genomics/microSALT.git`
* `conda config --add channels bioconda`
* `conda create -n microSALT python=3.6`
* `source activate microSALT`
* `conda install blast spades trimmomatic`
* `cd microSALT && pip install -r requirements.txt && pip install .`
* Perform all steps under section  __Configuration__

## Configuration
Copy the configuration file `configExample.json` to `~/.microSALT/config.json` _or_ place it wherever and point $MICROSALT_CONFIG to it.

Modify the line `SQLALCHEMY_DATABASE_URI` to fill out your database credentials. For production purposes, set the `DEBUG` flag to False.

Edit the other fields to match your environment.

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

### Genologics python3 bug fix
Change line 5 of `config.py` to `import configparser as ConfigParser` to fix the bug.
To find the path of the file, simply run `microSALT` and note where the log points to.

## Usage
* Use the `start` function to start sbatch job(s), producing output to `folders['results']`.
* After you have been informed of job completetion (through e-mail). Use the `finish` function to upload parsed results to the SQL back-end and produce reports (HTML & CSV).
* Various functionality, including adding new reference organisms and re-generating reports; are stored under the `util` command.
