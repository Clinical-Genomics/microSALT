<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

microbial Sequence Analysis and Loci-based Typing pipeline

The microbial sequence analysis and loci-based typing pipeline (microSALT) is used to analyse microbial samples.
It produces a quality control of the sample, determines a sample's organism specific sequence type, and its resistance pattern. microSALT also provides a database storage solution and report generation of these results.

microSALT uses a combination of python, sqLite and flask. Python is used for the majority of functionality, the database is handled through sqLite and the front-end is handled through flask. All analysis activity by microSALT requires a SLURM cluster.

## Requirements
### Hardware
* A slurm enabled HPC
* A (clarity) LIMS server
* A sqLite service

### Software
* Conda
* Python 3.6

## Quick installation
* `git clone https://github.com/Clinical-Genomics/microSALT.git`
* `cd microSALT && bash install.sh`
* Perform all steps under section  __Configuration__

## Configuration
Copy the configuration file `configExample.json` to `~/.microSALT/config.json` _or_ place it wherever and point $MICROSALT_CONFIG to it.

Edit the fields to match your environment.

### Genologics Configuration
_Genologics is likely already installed on your system. If such, this section can be skipped_
Create `$HOME/.genologicsrc` with the following formatting:
```
[genologics]
BASEURI=https://yourlims.corporation.se/
USERNAME=your_username
PASSWORD=your_password
[logging]
MAIN_LOG=/tmp/lims.log
```

#### Genologics python3 bug fix
Change line 5 of `config.py` to `import configparser as ConfigParser` to fix the bug.
To find the path of the file, simply run `microSALT` and note where the log points to.

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

## Credits
* Isak Sylvin - Lead developer
* Emma Sernstad - Accreditation ready reports
* Tanja Normark - Contamination analysis, various issues
* Maya Brandi - Various issues
