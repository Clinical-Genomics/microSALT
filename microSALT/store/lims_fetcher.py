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

  def __init__(self, log):
    self.data = {}
    self.lims = Lims(BASEURI, USERNAME, PASSWORD)
    self.logger = log

  def get_lims_project_info(self, cg_projid):
    project = Project(self.lims, id=cg_projid)
    try:
      self.data.update({'date_received': project.open_date,
                               'CG_ID_project': cg_projid,
                               'Customer_ID_project' : project.name})
    except KeyError as e:
      self.logger.warn("Unable to fetch LIMS info for project {}\nSource: {}".format(cg_projid, str(e)))

  def get_lims_sample_info(self, cg_sampleid):
    sample = Sample(self.lims, id=cg_sampleid)
    #TODO: Should control samples be analyzed?
    if 'Strain' not in sample.udf or sample.udf['Strain'] == 'Other':
      self.logger.warn("Unspecific strain specified for sample {}. Assuming control sample, thus ignoring."\
      .format(self.data['organism'], cg_sampleid))
    elif sample.udf['Strain'] == 'VRE':
      try:
        self.data.update({'organism' : sample.udf['Comment']})
      except Exception as e:
        self.logger.warn("Ambigious organism {} found in sample {}, but no comment clarifies.\nSource: {}"\
        .format(self.data['organism'], cg_sampleid, str(e)))
    else:
      try:
        self.data.update({'CG_ID_project': sample.project.id,
                             'CG_ID_sample': cg_sampleid,
                             'Customer_ID_sample' : sample.name,
                             'organism' : sample.udf['Strain']})
      except KeyError as e:
        self.logger.warn("Unable to fetch LIMS info for sample {}. Review LIMS data.\nSource: {}"\
        .format(cg_sampleid, str(e)))
