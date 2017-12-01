#!/usr/bin/env python

import click
import glob
import os
import re
import sys
import time
import yaml

from genologics.lims import Lims
# Should probably call these items directly since we're now up to 3 config files
from genologics.config import BASEURI,USERNAME,PASSWORD
from genologics.entities import Project, Sample

class LIMS_Fetcher():

  def __init__(self):
    self.data = {}
    self.lims = Lims(BASEURI, USERNAME, PASSWORD)

  def get_lims_project_info(self, cg_projid):
    project = Project(self.lims, id=cg_projid)
    try:
      self.data.update({'date_received': project.open_date,
                               'Customer_ID_project' : project.name})
    except KeyError:
      pass

  def get_lims_sample_info(self, cg_sampleid):
    sample = Sample(self.lims, id=cg_sampleid)
    try:
      self.data.update({'CG_ID_project': sample.project.id,
                             'Customer_ID_sample' : sample.name,
                             'organism' : sample.udf['Strain']})
    except KeyError:
      pass
