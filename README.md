[![Build status](https://github.com/clinical-genomics/microsalt/actions/workflows/run_tests.yml/badge.svg)](https://github.com/clinical-genomics/microsalt/actions/workflows/run_tests.yml)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.4026043-blue)](https://doi.org/10.5281/zenodo.4026043)

<p align="center">
  <a href="https://github.com/Clinical-Genomics/microSALT">
    <img width="1000" height="250" src="microSALT/artwork/microsalt.jpg"/>
  </a>
</p>

**Microbial Sequence Analysis and Loci-based Typing pipeline**

_The microbial sequence analysis and loci-based typing pipeline (microSALT) is
used to analyse microbial samples. It produces a quality control of the
sample, determines a sample's organism specific sequence type, and its
resistance pattern. microSALT also provides a database storage solution and
report generation of these results._

_microSALT uses a combination of Python, MySQL and Jinja2. Python is used for
the majority of functionality, the database is handled through MySQL via
SQLAlchemy and reports are rendered through Jinja2. All analysis activity by
microSALT requires a SLURM cluster._

## Installation

### Quick install

> [!IMPORTANT]
> This install requires `uv` to be installed on the system. For installation instructions, see [https://docs.astral.sh/uv/getting-started/installation/](https://docs.astral.sh/uv/getting-started/installation/).

`bash <(curl https://raw.githubusercontent.com/Clinical-Genomics/microSALT/master/install.sh)`

### Manual install

1. Clone the repository and enter the directory
1. Checkout the desired branch
1. install package using `uv pip install .`

## Configuration

Copy the configuration file anywhere and.

`cp configExample.json $HOME/.microSALT/config.json`

> [!IMPORTANT]
> **Then edit the fields to match your environment**.

## Installing containers

microSALT uses [Singularity](https://sylabs.io/singularity/) containers to run the various tools used in the analysis. These containers are available on Clinical Genomics' DockerHub, and can be pulled using the following command:

```Console
singularity pull docker://clinicalgenomics/microsalt-blast:latest
singularity pull docker://clinicalgenomics/microsalt-bwa:latest
singularity pull docker://clinicalgenomics/microsalt-picard:latest
singularity pull docker://clinicalgenomics/microsalt-quast:latest
singularity pull docker://clinicalgenomics/microsalt-samtools:latest
singularity pull docker://clinicalgenomics/microsalt-skesa:latest
singularity pull docker://clinicalgenomics/microsalt-trimmomatic:latest
```

> [!NOTE]
> Remember to enter the correct path to the singularity images in the configuration file.

## Usage

- `microsalt analyse` contains functions to start sbatch job(s) & produce
    output to `folders['results']`. Afterwards the parsed results are uploaded
    to the SQL back-end and produce reports (HTML), which are then automatically
    e-mailed to the user.
- `microsalt utils` contains various functionality, including generating the
    sample description JSON, manually adding new reference organisms and
    re-generating reports.

## Setup

Before running microSALT, the user must run the `setup` command, which will create the necessary database tables and download the necessary databases. This only needs to be run once, and can be run again if the user wants to reset the database or download new databases.

The setup is also dependent on

```Shell
microsalt setup
```

## Retrieving credentials

The credentials to access the [pubMLST and Pasteur](#mlst-definitions) database can be retrieved by running the following command:

```Shell
microsalt utils get_bigsdb_credentials
```

This will allow the user to specify which database they want to retrieve credentials for. Given that the user has given the correct information in the [Configuration section](#configuration), the credentials will be retrieved and stored on disk for later use.

## Databases

### MLST Definitions

microSALT will automatically download & use the MLST definitions for any
organism on [pubMLST](https://pubmlst.org/databases) or [Pasteur](https://bigsdb.pasteur.fr/). Other definitions may be
used, as long as they retain the same format.

### Resistance genes

microSALT will automatically download & use the resistance genes of [ResFinder](https://bitbucket.org/genomicepidemiology/resfinder).
Any definitions will work, as long as they retain the same formatting.

## Requirements

### Hardware

- A [SLURM](https://slurm.schedmd.com) enabled HPC system

### Software

- [uv](https://docs.astral.sh/uv) >= 0.4
- Python >= 3.10
- [MySQL](https://www.mysql.com) server

## Contributing to this repo

This repository follows the Github flow approach to adding updates.
For more information, see https://guides.github.com/introduction/flow/

## Credits

- Isak Sylvin - Lead developer
- Emma Sernstad - Accreditation ready reports
- Tanja Normark - Various issues
- Maya Brandi - Various issues
