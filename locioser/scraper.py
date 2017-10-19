"""This initial script scrapes output files for data and add them to the database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
import time
import yaml

import pdb # debug

with open("{}/config.yml".format(os.path.dirname(os.path.realpath(__file__))), 'r') as conf:
  config = yaml.load(conf)


#Take in BLAST file
#Take in assembly stats
#Take in QC yaml output (Robin knows)
#Take in LIMSid -> udfinfo (Kenny knows)

#Parse files
#Dump into database (via db manipulator)
