<p align="center">
  <a href="https://github.com/sylvinite/microSALT">
    <img width="1000" height="250" src="artwork/microsalt.jpg"/>
  </a>
</p>

# microSALT
microbial Sequence Analysis and Loci-based Typing pipeline

Currently a work in progress. 

First release estimated time of arrival by christmas.

## Dependencies
* Python 3.6
* MySQL

* BLAST
* CutAdapt
* SpaDES

* Genologics, `https://github.com/SciLifeLab/genologics.git`

## Configuration
### Flask/Database configuration
An additional config.py file is required to be present in the instance folder.

Formatting of the file is as follows:
```
# -*- coding: utf-8 -*-

SQLALCHEMY_DATABASE_URI= 'mysql+pymysql://DB_USER:DB_PASSWORD@DB_HOST:DB_PORT/DB_DATABASENAME'
SQLALCHEMY_TRACK_MODIFICATIONS= False
DEBUG= True
```
Debug statement has to be omitted for any production usage.
### Paths
Review `microSALT/config/paths_and_headers.yml` so that it accurately uses the headers and paths you desire.

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

Additionally replace `ConfigParser` with `configparser` in `config.py` for python3 support. 
