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

# TODO: Rewrite so samples use seperate objects
class Scraper():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.infolder = os.path.abspath(infolder)
    self.sampledir = ""
   
    last_folder = os.path.basename(os.path.normpath(self.infolder)) 
    self.name = last_folder.split('_')[0]
    #TODO: Replace date from dir with entry from analysis files/database
    self.date = "{} {}".format(re.sub('\.','-', last_folder.split('_')[1]), re.sub('\.',':', last_folder.split('_')[2]))
    self.db_pusher=DB_Manipulator(config, self.logger)
    self.lims_fetcher=LIMS_Fetcher(log, config)
    self.lims_sample_info = {}

  def scrape_project(self):
    #Scrape order matters a lot!
    self.lims_fetcher.get_lims_project_info(self.name)
    self.scrape_projectinfo()
    for dir in os.listdir(self.infolder):
     if os.path.isdir("{}/{}".format(self.infolder, dir)): 
       self.sampledir = "{}/{}".format(self.infolder, dir)
       self.name = dir
       self.lims_fetcher.get_lims_sample_info(dir)
       self.scrape_sampleinfo()
       self.scrape_all_loci()

  def scrape_sample(self):
    #Scrape order matters a lot!
    self.sampledir = self.infolder
    self.lims_fetcher.get_lims_sample_info(self.name)
    self.lims_fetcher.get_lims_project_info(self.lims_fetcher.data['CG_ID_project'])

    self.scrape_projectinfo()
    self.scrape_sampleinfo()
    self.scrape_all_loci()

  def scrape_all_loci(self):
    q_list = glob.glob("{}/loci_query_*".format(self.sampledir))
    for file in q_list:
      self.scrape_single_loci(file)
    #Requires all loci results to be initialized
    try:
      ST = self.db_pusher.alleles2st(self.name)
      self.db_pusher.upd_rec({'CG_ID_sample':self.name}, 'Samples', {'ST':ST})
    except Exception as e:
      pass

  def scrape_projectinfo(self):
    proj_col=dict()
    proj_col['CG_ID_project'] = self.name
    proj_col['Customer_ID_project'] = self.lims_fetcher.data['Customer_ID_project']
    proj_col['date_ordered'] = self.lims_fetcher.data['date_received']
    self.db_pusher.add_rec(proj_col, 'Projects')

  def scrape_sampleinfo(self):
    """Identifies sample values"""
    try:
      sample_col = self.db_pusher.get_columns('Samples') 
      sample_col['CG_ID_sample'] = self.lims_fetcher.data['CG_ID_sample']
      sample_col['CG_ID_project'] = self.lims_fetcher.data['CG_ID_project']
      sample_col['Customer_ID_sample'] = self.lims_fetcher.data['Customer_ID_sample']
      sample_col["date_analysis"] = self.date
      self.db_pusher.add_rec(sample_col, 'Samples')
    except Exception as e:
      #Ignores unloaded samples
      pass

  def scrape_single_loci(self, infile):
    """Scrapes a single blast output file for MLST results"""
    seq_col = self.db_pusher.get_columns('Seq_types') 
    if not os.path.exists(self.sampledir):
      self.logger.error("Invalid file path to infolder, {}".format(self.sampledir))
      sys.exit()
    try:
      with open("{}".format(infile), 'r') as insample:
        organism = self.lims_fetcher.get_organism_refname(self.name)
        self.db_pusher.upd_rec({'CG_ID_sample' : self.name}, 'Samples', {'organism': organism})
        seq_col["CG_ID_sample"] = self.name

        for line in insample:
          #Ignore commented fields
          if not line[0] == '#':
            elem_list = line.rstrip().split("\t")
            if not elem_list[1] == 'N/A':
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
    except Exception as e:
      pass
