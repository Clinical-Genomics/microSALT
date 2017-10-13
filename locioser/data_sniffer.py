"""Initial script to fetch data from QC, Assembly and BLAST output and put them into a database.
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


