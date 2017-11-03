"""This initial script scrapes output files for data and add them to the database
   Heavy WIP
   By: Isak Sylvin, @sylvinite"""

#!/usr/bin/env python

import click
import os
import re
import time
import yaml

from microSALT import db_manipulator

import pdb # debug

class Scraper():

  def __init__(self, infile, config):
    self.config = config
    self.infile = os.path.abspath(infile)
    self.blastdb=db_manipulator.DB_Manipulator()

  def scrape_blast_loci(self):
    #Assign each element correct keys
    self.blastdb.create_blasttable()
    columns = self.blastdb.get_blastcolumns() 

    with open("{}".format(self.infile), 'r') as inblast:
      for line in inblast:
        #Ignore commented fields
        if not line[0] == '#':
          elem_list = line.rstrip().split("\t")
          # remove " + " elem 0
          elem_list[0] = re.search('(\w+)', elem_list[0]).group(0)
          # split elem 2 into contig node_NO, length, cov
          nodeinfo = elem_list.pop(2).split('_')
          elem_list.insert(2, "{}_{}".format(nodeinfo[0], nodeinfo[1]))
          elem_list.insert(3, nodeinfo[3])
          elem_list.insert(4, nodeinfo[5]) 
        
          #Get run and date from folder structure
          rundir = os.path.basename(os.path.dirname(self.infile))
          rundir = rundir.split('_')
          elem_list.insert(0, rundir[0])
          rundir[1] = re.sub('\.','-',rundir[1])
          rundir[2] = re.sub('\.',':',rundir[2])
          elem_list.insert(1, "{} {}".format(rundir[1], rundir[2]))

          #Since everything is ordered correctly, match values to keys
          i=0
          inrecord = dict()
          while i < len(elem_list):
            inrecord[columns[i]] = elem_list[i]
            i = i+1
          #Needs to check if entry already exists
          #pdb.set_trace()
          self.blastdb.add_blastrecord(inrecord)

#Take in assembly stats
#Take in QC yaml output (Robin knows)
#Take in LIMSid -> udfinfo (Kenny knows)

#Parse files
#Dump into database (via db manipulator)

