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
* `conda config --add channels bioconda`
* `conda config --add channels conda-forge`
* `conda create -n MLST python=3.6`
* `conda install blast trimmomatic spades`
* `source activate MLST`
* `git clone https://github.com/sylvinite/microSALT.git`
* Perform all steps under section  __Configuration__
* `cd microSALT && pip install -r requirements.txt && pip install -e . && cd ..`

## Configuration
Rename the configuration file under folder `instance` to `config.json` and modify the line `SQLALCHEMY_DATABASE_URI` to correctly point to your database. For production purposes, set the `DEBUG` flag to False. Review the paths in the file to accurately represent the file paths and sbatch headers required by your system.

### LIMS
Create `$HOME/.genologicsrc` with the following formatting:
```
[genologics]
BASEURI=https://yourlims.example.com:8443
USERNAME=your_username
PASSWORD=your_password
[logging]
MAIN_LOG=/home/glsai/your_main_log_file
```

### Genologics python3 support
Change line 5 of `config.py` to `import configparser as ConfigParser` for python3 support.
Easiest way to find your `config.py` file is to run `microSALT` after all other configuration steps have been executed.

## Usage
* Use the `create` function to generate sbatch job(s) defined under `folders['results']`. Manually start them via the `concatinated.sh` script.
* After the jobs have been finished. Use the `scrape` function to upload parsed results to the SQL backend.
* Use the `view` function to start a flask instance to view the results. Point your browser at `http://127.0.0.1:5000/microSALT`
* Navigate to your run, and print the results to PDF format (Command/Ctrl + P) if requested.
