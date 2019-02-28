#!/usr/bin/env python

import os
import re

from datetime import datetime
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
    samplelist = self.samples_in_project(cg_projid)

    custids = list()
    for sample in samplelist:
      self.load_lims_sample_info(sample)
      custids.append(self.data['Customer_ID'])
    if not custids[1:] == custids[:-1]:
      raise Exception("Project {} contains multiple Customer IDs".format(cg_projid))
    try:
      #Resolves old format
      if ' ' in project.name:
        realname = project.name.split(' ')[0]
      else:
        realname = project.name

      tmp = project.open_date.split('-')
      newdate = datetime(int(tmp[0]), int(tmp[1]), int(tmp[2]))
      self.data.update({'date_received': newdate,
                               'CG_ID_project': cg_projid,
                               'Customer_ID_project' : realname,
                               'Customer_ID': custids[0]})
    except KeyError as e:
      self.logger.warn("Unable to fetch LIMS info for project {}\nSource: {}".format(cg_projid, str(e)))

  def samples_in_project(self, cg_projid, external=False):
    """ Returns a list of sample names for a project"""
    output = list()
    if not external:
      samples = self.lims.get_samples(projectlimsid=cg_projid)
    else:
      samples = self.lims.get_samples(projectname=cg_projid)
    for s in samples:
      output.append(s.id)
    return output

  def load_lims_sample_info(self, cg_sampleid, external=False):
    """ Loads all utilized LIMS info. Organism assumed to be written as binomial name """
    libprep_date = ""
    seq_date = ""
    if external:
      sample = self.lims.get_samples(name=cg_sampleid)
      if len(sample) != 1:
        self.logger.error("Sample ID {} resolves to multiple entries".format(cg_sampleid))
      sample = sample[0]
      seq_date = self.get_date(sample.id,type="sequencing")
      libprep_date = self.get_date(sample.id,type="libprep")
    else:
      try:
        sample = Sample(self.lims, id=cg_sampleid)
        seq_date = self.get_date(cg_sampleid,type="sequencing")
        libprep_date = self.get_date(cg_sampleid,type="libprep")
      except Exception as e:
        self.logger.error("LIMS connection timeout")
    organism = "Unset"
    if 'Strain' in sample.udf and organism == "Unset":
      #Backwards compat, MUST hit first
      if sample.udf['Strain'] == 'VRE':
        if 'Reference Genome Microbial' in sample.udf:
          if sample.udf['Reference Genome Microbial'] == 'NC_017960.1':
            organism = 'Enterococcus faecium'
          elif sample.udf['Reference Genome Microbial'] == 'NC_004668.1':
            organism = 'Enterococcus faecalis'
        elif 'Comment' in sample.udf:
          organism = sample.udf['Comment']
      elif sample.udf['Strain'] != 'Other' and sample.udf['Strain'] != 'other':
        organism = sample.udf['Strain']
      elif (sample.udf['Strain'] == 'Other' or sample.udf['Strain'] == 'other') and 'Other species' in sample.udf:
        organism = sample.udf['Other species']
    if 'Reference Genome Microbial' in sample.udf and organism == "Unset":
      if sample.udf['Reference Genome Microbial'] == 'NC_002163':
        organism = "Campylobacter jejuni"
      elif sample.udf['Reference Genome Microbial'] == 'NZ_CP007557.1':
        organism = 'Klebsiella oxytoca'
      elif sample.udf['Reference Genome Microbial'] == 'NC_000913.3':
        organism = 'Citrobacter freundii'
      elif sample.udf['Reference Genome Microbial'] == 'NC_002516.2':
        organism = 'Pseudomonas aeruginosa'
    elif 'Comment' in sample.udf and organism == "Unset":
      organism = sample.udf['Comment']
    # Consistent safe-guard
    elif organism == "Unset":
      organism = "Other"
      self.logger.warn("Unable to resolve ambigious organism found in sample {}."\
      .format(cg_sampleid))
    try:
      self.data.update({'CG_ID_project': sample.project.id,
                           'CG_ID_sample': sample.id,
                           'Customer_ID_sample' : sample.name,
                           'organism' : organism,
                           'priority' : sample.udf['priority'],
                           'Customer_ID': sample.udf['customer'],
                           'application_tag': sample.udf['Sequencing Analysis'],
                           'date_sequencing': seq_date,
                           'date_libprep': libprep_date
})
    except KeyError as e:
      self.logger.warn("Unable to fetch LIMS info for sample {}. Review LIMS data.\nSource: {}"\
      .format(cg_sampleid, str(e)))

  def get_organism_refname(self, sample_name, external=False):
    """Finds which reference contains the same words as the LIMS reference
       and returns it in a format for database calls."""
    self.load_lims_sample_info(sample_name, external)
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
      self.logger.warn("Unable to find existing reference for {}, strain {} has no reference match\nSource: {}"\
      .format(sample_name, lims_organ, e))

#  def get_date(self, sample_id, type=""):
#    """ Returns the most recent sequencing date of a sample """
#    if type == "sequencing":
#        steps = ["CG002 - Illumina Sequencing (Illumina SBS)", "CG002 Illumina SBS (HiSeq X)"]
#    elif type == "libprep":
#        steps = "CG002 - Aggregate QC (Library Validation)"
#    else:
#        raise Exception("Attempted to get date for {} but no step defined".format(sample_id))
#    try:
#        arts = self.lims.get_artifacts(samplelimsid = sample_id, process_type = steps)
#        date_list = [a.parent_process.date_run for a in arts]
#    except Exception as e:
#        pass
#    return max(date_list)

  def get_date(self, sample_id, type=""):
    """ Returns the most recent sequencing date of a sample """
    date_list = list()
    if type == "sequencing":
      steps = ["CG002 - Illumina Sequencing (Illumina SBS)", "CG002 Illumina SBS (HiSeq X)"]
    elif type == "libprep":
      steps = ["CG002 - Aggregate QC (Library Validation)"]
    else:
      raise Exception("Attempted to get date for {} but no step defined".format(sample_id))
    for step in steps:
      try:
        arts = self.lims.get_artifacts(samplelimsid = sample_id, process_type = step)
        date_list = date_list + [a.parent_process.date_run for a in arts]
      except Exception as e:
        pass
    dp = max(date_list).split('-')
    return datetime(int(dp[0]), int(dp[1]), int(dp[2]))

#DEBUG: Implement later!
#   def get_method_document(sample_id, first_date = True):
#     steps = ['CG002 - Cluster Generation (HiSeq X)', 'CG002 - Cluster Generation (Illumina SBS)']
#
#     arts = lims.get_artifacts(samplelimsid = sample_id, process_type = steps)
#     processes = [(a.parent_process.date_run, a.parent_process.date_run) for a in arts]
#     processes = list(set(processes))
#     if first_date
#       process = sorted(processes)[0]
#     else:
#       process = sorted(processes)[-1]
#     dokument = process.udf.get('Method Document')
#     version = process.udf.get('document version')
#     return dokument, version
