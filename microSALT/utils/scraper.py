"""This initial script scrapes output files for data and add them to the database
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

class Scraper():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.infolder = os.path.abspath(infolder)
    
    self.db_pusher=DB_Manipulator(config, self.logger)
    self.lims_fetcher=LIMS_Fetcher()
    self.CG_ID_sample = ""
    self.sampledir = ""
    self.lims_sample_info = {}

  def scrape_project(self):
    #Scrape order matters a lot!
    self.lims_fetcher.get_lims_project_info(os.path.basename(os.path.normpath(infolder)))
    self.scrape_projectinfo()
    for (dirpath, dirnames, filenames) in os.walk(indir):
      for dir in dirnames:
        self.scrape_sampleinfo()
        self.scrape_all_loci()

  def scrape_sample(self):
    #Scrape order matters a lot!
    self.sampledir = os.path.basename(os.path.normpath(self.infolder))
    self.CG_ID_sample = self.sampledir.split('_')[0]
    self.lims_fetcher.get_lims_sample_info(self.CG_ID_sample)
    self.lims_fetcher.get_lims_project_info(self.lims_fetcher.data['CG_ID_project'])

    self.scrape_projectinfo()
    self.scrape_sampleinfo()
    self.scrape_all_loci()

  def scrape_all_loci(self):
    q_list = glob.glob("{}/loci_query_*".format(self.infolder))
    for file in q_list:
      self.scrape_single_loci(file)
    #Requires all loci results to be initialized
    ST = self.db_pusher.alleles2st(self.CG_ID_sample)
    self.db_pusher.upd_rec_orm({'CG_ID_sample':self.CG_ID_sample}, 'Samples', {'ST':ST})

  def scrape_projectinfo(self):
    proj_col=dict()
    proj_col['CG_ID_project'] = self.lims_fetcher.data['CG_ID_project']
    proj_col['Customer_ID_project'] = self.lims_fetcher.data['Customer_ID_project']
    proj_col['date_ordered'] = self.lims_fetcher.data['date_received']
    self.db_pusher.add_rec_orm(proj_col, 'Projects')

  def scrape_sampleinfo(self):
    """Identifies sample values"""
    sample_col = self.db_pusher.get_columns_orm('Samples')
    sample_col["CG_ID_sample"] = self.CG_ID_sample
    rundir = self.sampledir.split('_')
    rundir[1] = re.sub('\.','-',rundir[1])
    rundir[2] = re.sub('\.',':',rundir[2])

    sample_col['CG_ID_project'] = self.lims_fetcher.data['CG_ID_project']
    sample_col['Customer_ID_sample'] = self.lims_fetcher.data['Customer_ID_sample']
    sample_col["date_analysis"] = "{} {}".format(rundir[1], rundir[2])
    self.db_pusher.add_rec_orm(sample_col, 'Samples')

  #TODO: Look over column assignment, since new input probably screwed things over
  def scrape_single_loci(self, infile):
    """Scrapes a single blast output file for MLST results"""
    seq_col = self.db_pusher.get_columns_orm('Seq_types') 
    if not os.path.exists(self.infolder):
      self.logger.error("Invalid file path to infolder, {}".format(self.infolder))
      sys.exit()
    with open("{}".format(infile), 'r') as insample:
      insample.readline()
      insample.readline()
      #Takes ref db/organism from resultfile. Awkward.
      db = insample.readline().rstrip().split(' ')
      organ = db[2].split('/')[-2]
      self.db_pusher.upd_rec_orm({'CG_ID_sample' : self.CG_ID_sample}, 'Samples', {'organism': organ})

      seq_col["CG_ID_sample"] = self.CG_ID_sample

      for line in insample:
        #Ignore commented fields
        if not line[0] == '#':
          elem_list = line.rstrip().split("\t")

          seq_col["identity"] = elem_list[4]
          seq_col["evalue"] = elem_list[5]
          seq_col["bitscore"] = elem_list[6]
          seq_col["contig_start"] = elem_list[7]
          seq_col["contig_end"] = elem_list[8]
          seq_col["loci_start"] = elem_list[9]
          seq_col["loci_end"] =  elem_list[10]
          seq_col["haplotype"] = elem_list[1]

         
          # Split elem 3 in loci (name) and allele (number) 
          seq_col["loci"] = elem_list[3].split('_')[0]
          seq_col["allele"] = int(elem_list[3].split('_')[1])

          # split elem 2 into contig node_NO, length, cov
          nodeinfo = elem_list[2].split('_')
          seq_col["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
          seq_col["contig_length"] = nodeinfo[3]
          seq_col["contig_coverage"] = nodeinfo[5]

          self.db_pusher.add_rec_orm(seq_col, 'Seq_types')


    self.logger.info("Added a record to the database")
