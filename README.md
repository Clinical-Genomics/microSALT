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

### QC-only analysis

To run read/alignment QC without assembly, MLST, resistance or ExPEC analysis:

```Shell
microsalt --config /path/to/config.json analyse --qc-only --input /path/to/fastq sampleinfo.json
```

`--qc_only` is also accepted. This applies to both single-sample and project runs.
The final job automatically uses `utils finish --track qc --report qc`, skipping
typing/resistance scraping and generating only the QC report. A configured
reference genome is still required for alignment. Without this flag, the full
analysis and reports are unchanged.

### Analysis validation checks

Before submitting any SLURM jobs, `analyse` checks sample metadata, paired gzip
FASTQ inputs, reference genome FASTA files and indexes, and (for full analysis)
MLST and resistance resources. It then prepares every sample script and the
reporting inputs. Setup errors are collected across samples and references and
reported together; any setup error exits nonzero without submitting sample or
reporting jobs. `--dry` uses the same checks but never submits jobs.

Whitespace is removed from sample-info `reference` values before lookup,
downloading and job generation. The original sample-info file is not modified;
the run's saved sample info contains the normalized reference values.

Downloaded genomes must contain valid, non-empty nucleotide FASTA records and
index successfully before they are published. Existing cached FASTA files are
also validated, so an older NCBI error response saved as `.fasta` blocks the run
with an actionable error. Inspect and repair that file before retrying.
NCBI assembly accessions (`GCF_...` / `GCA_...`) are looked up in the Assembly
database. microSALT downloads the complete genomic FASTA (including all contigs
and plasmids) from the matching RefSeq or GenBank assembly directory over HTTPS.
Versioned accessions must match the requested version; microSALT does not silently
substitute a newer version or the paired assembly. Nucleotide accessions still use
the nucleotide database. Both download paths use the same validation and indexing
checks, and keep the supplied accession as the cached `<reference>.fasta` filename.
You can also install a valid local `<reference>.fasta` in `folders.genomes`.

`--skip_update` disables reference downloading and indexing, **not validation**:
all required reference files and current indexes must already exist. QC-only
runs still require an alignment genome, but do not require typing resources.

An organism without a published MLST database or scheme is not a setup error.
microSALT warns and skips MLST searches and scraping for that organism, while
retaining assembly, resistance analysis and QC. Missing MLST version dates only
disable version comparison; available schemes can still be downloaded. Network,
download, indexing and corrupt cached-resource failures still block submission.
With `--skip_update`, an organism with no cached scheme is likewise analysed
without MLST rather than requiring an online availability check.

If SLURM rejects a submission after setup succeeds, the command reports the
failure and any already accepted analysis job IDs. This is distinct from a
setup failure: those accepted jobs may already be running.

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
