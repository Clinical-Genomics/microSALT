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

  def __init__(self, config, log):
    self.data = {}
    self.lims = Lims(BASEURI, USERNAME, PASSWORD)
    self.logger = log
    self.config = config

  def load_lims_project_info(self, cg_projid):
    project = Project(self.lims, id=cg_projid)
    try:
      self.data.update({'date_received': project.open_date,
                               'CG_ID_project': cg_projid,
                               'Customer_ID_project' : project.name})
    except KeyError as e:
      self.logger.warn("Unable to fetch LIMS info for project {}\nSource: {}".format(cg_projid, str(e)))

  def load_lims_sample_info(self, cg_sampleid):
    sample = Sample(self.lims, id=cg_sampleid)
    if 'Strain' in sample.udf: 
      organism = sample.udf['Strain']
      if sample.udf['Strain'] == 'VRE':
        if 'Reference Genome Microbial' in sample.udf:
          if sample.udf['Reference Genome Microbial'] == 'NC_017960.1':
            organism = 'Enterococcus faecium'
          elif sample.udf['Reference Genome Microbial'] == 'NC_004668.1':
            organism = 'Enterococcus faecalis'
        elif 'Comment' in sample.udf:
          organism = sample.udf['Comment']
    else:
      #Consistency
      organism = "other"
      self.logger.warn("Unable to resolve ambigious organism found in sample {}."\
      .format(cg_sampleid))
    try:
      self.data.update({'CG_ID_project': sample.project.id,
                           'CG_ID_sample': cg_sampleid,
                           'Customer_ID_sample' : sample.name,
                           'organism' : organism})
    except KeyError as e:
      self.logger.warn("Unable to fetch LIMS info for sample {}. Review LIMS data.\nSource: {}"\
      .format(cg_sampleid, str(e)))

  def get_organism_refname(self, sample_name):
    """Gets organism/strain name, returns formatted name for mlst reference use"""
    self.load_lims_sample_info(sample_name)
    lims_organ = self.data['organism'].lower()
    orgs = os.listdir(self.config["folders"]["references"])
    organism = re.split('\W+', lims_organ)
    try:
      refs = 0
      for target in orgs:
        hit = 0
        for piece in organism:
          if piece in target:
            hit +=1
          #For when people misspell the strain in the orderform
          elif piece == "pneumonsiae" and "pneumoniae" in target:
            hit +=1
          else:
            break
        if hit == len(organism):
          return target
    except Exception as e:
      self.logger.warn("Unable to find reference for {}, strain {} has no reference match\nSource: {}".format(sample_name, lims_organ, e))

