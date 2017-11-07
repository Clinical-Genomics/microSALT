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
    self.blastdb=db_manipulator.DB_Manipulator(config)

  def scrape_blast_loci(self):
    #Assign each element correct keys
    columns = self.blastdb.get_blastcolumns() 


    with open("{}".format(self.infile), 'r') as inblast:
      #Grab database = organism from blast output
      inblast.readline()
      inblast.readline()
      db = inblast.readline()
      db = db.rstrip().split(' ')
      columns['organism'] = os.path.basename(os.path.normpath(db[2]))
      
      for line in inblast:
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

          #Get run and date from folder structure
          rundir = os.path.basename(os.path.dirname(self.infile))
          rundir = rundir.split('_')
          columns["run"] = rundir[0]
          rundir[1] = re.sub('\.','-',rundir[1])
          rundir[2] = re.sub('\.',':',rundir[2])
          columns["date_analysis"] = "{} {}".format(rundir[1], rundir[2])

          #Get allele from ST number
          columns['allele'] = self.blastdb.st2allele(columns) 
          self.blastdb.add_blastrecord(columns)

#Take in assembly stats
#Take in QC yaml output (Robin knows)
#Take in LIMSid -> udfinfo (Kenny knows)

#Parse files
#Dump into database (via db manipulator)

