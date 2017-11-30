#!/usr/bin/env python

import click
import glob
import os
import pdb
import re
import sys
import time
import yaml

from genologics.lims import Lims
# Should probably call these items directly since we're now up to 3 config files
from genologics.config import BASEURI,USERNAME,PASSWORD
from genologics.entities import Sample

class LIMS_Fetcher():

  def __init__(self):
    self.data = {}

  def get_lims_sample_info(self, cg_sampleid):
    lims = Lims(BASEURI, USERNAME, PASSWORD)
    sample = Sample(lims, id=cg_sampleid)
    self.data = {'date_completed' : sample.date_completed,
                             'date_received' : sample.date_received,
                             'CG_ID_project': sample.project.id,
                             'Customer_ID_sample' : sample.name,
                             'Customer_ID_project' : sample.project.name}

