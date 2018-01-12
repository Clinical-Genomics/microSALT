"""Compares existing organism references with available and updates as needed
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import glob
import os
import re
import sys
import time
import yaml

from microSALT.store.db_manipulator import DB_Manipulator
from microSALT.store.lims_fetcher import LIMS_Fetcher

class Ref_Updater():

  def __init__(self, config, log):
    self.config = config
    self.logger = log

  def update_refs():
    pass
# Connect to resource

# Check version, names of exisitng references

# Check version of remote references

# Replace as needed


