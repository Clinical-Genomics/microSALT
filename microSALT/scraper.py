"""This initial script scrapes output files for data and add them to the database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import glob
import os
import pdb
import re
import sys
import time
import yaml

from microSALT import db_manipulator


class Scraper():

  def __init__(self, infolder, config, log):
    self.config = config
    self.logger = log
    self.infolder = os.path.abspath(infolder)
    self.db_pusher=db_manipulator.DB_Manipulator(config, self.logger)
    self.CG_ID_sample = ""

  def scrape_all_loci(self):
    q_list = glob.glob("{}/loci_query_*".format(self.infolder))
    #Due to FK samples MUST be scraped before seq_types
    self.scrape_sampleinfo()
    for file in q_list:
      #Unoptimized for multi. Calls sample once per seq_type
      self.scrape_single_loci(file)
    #Requires all loci results to be initialized
    ST = self.db_pusher.alleles2st(self.CG_ID_sample)
    self.db_pusher.upd_rec_orm({'CG_ID_sample':self.CG_ID_sample}, 'Samples', {'ST':ST})

  def scrape_sampleinfo(self):
    """Identifies sample values"""
    pcolumns = self.db_pusher.get_columns_orm('Samples')

    rundir = os.path.basename(os.path.normpath(self.infolder)).split('_')
    self.CG_ID_sample = rundir[0]
    pcolumns["CG_ID_sample"] = self.CG_ID_sample
    rundir[1] = re.sub('\.','-',rundir[1])
    rundir[2] = re.sub('\.',':',rundir[2])

    pcolumns['CG_ID_project'] = "P-{}".format(pcolumns["CG_ID_sample"])
    pcolumns["date_analysis"] = "{} {}".format(rundir[1], rundir[2])
    self.db_pusher.add_rec_orm(pcolumns, 'Samples')

  #TODO: Look over column assignment, since new input probably screwed things over
  def scrape_single_loci(self, infile):
    """Scrapes a single blast output file for MLST results"""
    columns = self.db_pusher.get_columns_orm('Seq_types') 
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

      columns["CG_ID_sample"] = self.CG_ID_sample

      for line in insample:
        #Ignore commented fields
        if not line[0] == '#':
          elem_list = line.rstrip().split("\t")

          columns["identity"] = elem_list[4]
          columns["evalue"] = elem_list[5]
          columns["bitscore"] = elem_list[6]
          columns["contig_start"] = elem_list[7]
          columns["contig_end"] = elem_list[8]
          columns["loci_start"] = elem_list[9]
          columns["loci_end"] =  elem_list[10]
          columns["haplotype"] = elem_list[1]

         
          # Split elem 3 in loci (name) and allele (number) 
          columns["loci"] = elem_list[3].split('_')[0]
          columns["allele"] = int(elem_list[3].split('_')[1])

          # split elem 2 into contig node_NO, length, cov
          nodeinfo = elem_list[2].split('_')
          columns["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
          columns["contig_length"] = nodeinfo[3]
          columns["contig_coverage"] = nodeinfo[5]

          self.db_pusher.add_rec_orm(columns, 'Seq_types')


    self.logger.info("Added a record to the database")
