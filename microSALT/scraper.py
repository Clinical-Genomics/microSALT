"""This initial script scrapes output files for data and add them to the database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import pdb
import re
import sys
import time
import yaml

from microSALT import db_manipulator


class Scraper():

  def __init__(self, infile, config, log):
    self.config = config
    self.logger = log
    self.infile = os.path.abspath(infile)
    self.db_pusher=db_manipulator.DB_Manipulator(config, self.logger)

  def scrape_loci_output(self):
    #Assign each element correct keys
    columns = self.db_pusher.get_columns_orm('Seq_types') 
    pcolumns = self.db_pusher.get_columns_orm('Samples')
    if not os.path.exists(self.infile):
      self.logger.error("Invalid file path to infile, {}".format(self.infile))
      sys.exit()
    with open("{}".format(self.infile), 'r') as insample:
      #Grab database = organism from sample output
      insample.readline()
      insample.readline()
      db = insample.readline()
      db = db.rstrip().split(' ')

      #Get run and date from folder structure
      rundir = os.path.basename(os.path.dirname(self.infile))
      rundir = rundir.split('_')
      pcolumns["CG_ID_sample"] = rundir[0]
      columns["CG_ID_sample"] = rundir[0]
      rundir[1] = re.sub('\.','-',rundir[1])
      rundir[2] = re.sub('\.',':',rundir[2])
 
      #TODO: Fetch from LIMS, setting placeholder for now
      pcolumns["CG_ID_project"] = "P-{}".format(columns["CG_ID_sample"])
      pcolumns["date_analysis"] = "{} {}".format(rundir[1], rundir[2])
      pcolumns['organism'] = os.path.basename(os.path.normpath(db[2]))
      pcolumns['ST'] = -1 #Should not be necessary, but safer.
      self.db_pusher.add_rec_orm(pcolumns, 'Samples')     

      for line in insample:
        #Ignore commented fields
        if not line[0] == '#':
          elem_list = line.rstrip().split("\t")

          columns["identity"] = elem_list[4]
          columns["evalue"] = elem_list[5]
          columns["bitscore"] = int(elem_list[6])
          columns["contig_start"] = elem_list[7]
          columns["contig_end"] = elem_list[8]
          columns["loci_start"] = elem_list[9]
          columns["loci_end"] =  elem_list[10]
          columns["haplotype"] = elem_list[1]

          #Remove the bp index from the ST hit
          #Note this ST is assumed, since one allele appear in many ST
          columns["assumed_ST"] = int(re.search('\w+', elem_list[3]).group(0)) 

          # remove " + " elem 0
          elem_list[0] = re.search('(\w+)', elem_list[0]).group(0)
          columns["loci"] = elem_list[0]
          # split elem 2 into contig node_NO, length, cov
          nodeinfo = elem_list.pop(2).split('_')
          columns["contig_name"] = "{}_{}".format(nodeinfo[0], nodeinfo[1])
          columns["contig_length"] = nodeinfo[3]
          columns["contig_coverage"] = nodeinfo[5]

          #TODO: REFIX THIS
          #Get allele from ST number
          columns['allele'] = self.db_pusher.st2allele(pcolumns['organism'], columns['loci'], columns['assumed_ST'])
          self.db_pusher.add_rec_orm(columns, 'Seq_types')

      #TODO: Fetch from LIMS, setting placeholder for now
      ST = self.db_pusher.alleles2st(columns['CG_ID_sample'])
      self.db_pusher.upd_rec_orm(pcolumns, 'Samples', {'ST':ST})
    self.logger.info("Added a record to the database")

