<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

microbial Sequence Analysis and Loci-based Typing pipeline
## Requirements
### Hardware
* A slurm enabled HPC
* A mySQL server

###Software
* Conda
* Python 3.6

## Quick installation
### Conda depdenency resolution
`conda config --add channels bioconda`
`conda config --add channels conda-forge`
`conda create -n MLST python=3.6`
`conda install blast trimmomatic spades`
`source activate MLST`

`git clone https://github.com/sylvinite/microSALT.git`
`git clone https://github.com/SciLifeLab/genologics.git`
Perform all steps under __Configuration__
`cd genologics && pip install -r requirements.txt && pip install . && cd ..`
`cd microSALT && pip install -r requirements.txt && pip install -e . && cd ..`

## Configuration
### Flask/Database configuration
Rename the configuration file under folder `instance` to `sqlalchemy_config.py` and modify the line `SQLALCHEMY_DATABASE_URI` to correctly point to your database. For production purposes, set the `DEBUG` flag to False.

### Paths file
Rename the configuration file under folder `instance` to `paths.yml`. Review the file to accurately represent the file paths and sbatch headers required by your system.

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
`ConfigParser` in `config.py` needs to be replaced with `configparser` for python3 support.

## Usage
* Use the `create` function to generate sbatch job(s) defined under `folders['results']`. Manually start them via the `concatinated.sh` script.
* After the jobs have been finished. Use the `scrape` function to upload parsed results to the SQL backend.
* Use the `view` function to start a flask instance to view the results. Point your browser at `http://127.0.0.1:5000/microSALT`
* Navigate to your run, and print the results to PDF format (Command/Ctrl + P) if requested.
